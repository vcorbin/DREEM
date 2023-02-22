#!/usr/bin/python

import numpy as np
import plotly.graph_objs as go
from plotly import tools
import pandas as pd
import logging
import plotly.offline as oplot
import EM_classes as EMClass

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import seaborn as sns
log = logging.getLogger(__name__)

class clusterPlots(object):
    def __init__(self, cluster_group, ref=None, sequence=None):
        """
        Creates a new instance of a clusterPlots
        :param cluster_group: Object of class EM_classes.clusterGroup
        """
        print('\n')

        if sequence is None:
            assert len(cluster_group.cluster_info.refs.keys()) == 1 or ref in cluster_group.cluster_info.refs.keys(), \
                "There are more than 1 reference in the cluster file comments section. You either have not specified " \
                "which reference you wish to use as a sequence or the reference you specified is not present in the " \
                "file."
            if ref is not None:
                if len(cluster_group.cluster_info.refs.keys()) == 1 and \
                                ref not in cluster_group.cluster_info.refs.keys():
                    print("Warning: there is only one reference listed in the cluster file (%s) and it does not " \
                          "correspond to the one you specified (%s). The former will be used instead.\n" \
                          % (cluster_group.cluster_info.refs.keys()[0], ref))
                    ref = cluster_group.cluster_info.refs.keys()[0]
                assert ref in cluster_group.cluster_info.refs.keys(), "The reference you specified (%s) is not " \
                                                                      "present in the comments section of the cluster" \
                                                                      " file." % ref
            else:
                assert len(cluster_group.cluster_info.refs.keys()) == 1, "There are more than 1 reference in the " \
                                                                         "cluster file comments section (%s) but you " \
                                                                         "have not specified which one you want to " \
                                                                         "use to plot the reference sequence with " \
                                                                         "the option -R." % \
                                                                    ", ".join(cluster_group.cluster_info.refs.keys())
                ref = list(cluster_group.cluster_info.refs.keys())[0]
            sequence = cluster_group.cluster_info.refs[ref]['seq']
        else:
            assert len(sequence) == cluster_group.df.shape[0], "Error: The sequence specified doesn't have the length" \
                                                               " (%s vs %s). Don't specify sequence in order to use" \
                                                               " the same one as the clustering." \
                                                               % (len(sequence), cluster_group.df.shape[0])

        print("sequence:", sequence)
        print('')

        cluster_df = cluster_group.df

        #assert len(sequence) == len(cluster_group.df), \
        #    "ERROR: The reference sequence doesn't have the same length as the clusters (%d vs %d)." % \
        #    (len(sequence), len(cluster_group.df))

        print(len(cluster_df), "length of cluster")
        print(len(sequence), "length of sequence")
        cluster_df['Base'] = np.array(list(sequence))
        cluster_df['Position'] = cluster_df.index

        print(cluster_df.head())

        self.num_clusters = cluster_df.shape[1] - 2
        self.clusters = cluster_df
        self.name = 'Clustering'
        self.filename = cluster_group.filename

    
    def setName(self, name):
        self.name = name
    
    def getNumClusters(self):
        return self.num_clusters
    
    def getClusters(self):
        return self.clusters
        
    def getBIC(self):
        pass # TODO:

    def saveResult(self):
        return self.clusters.to_csv(sep='\t', index=False)
        
    def linePlot(self, title=None):
        if title is None:
            title = 'Subpopulations within Region on '+self.name
        data = []
        for j in range(self.num_clusters):
            trace = go.Scatter(
                x=self.clusters['Position'],
                y=self.clusters['Cluster_'+str(j)],
                text=self.clusters['Base'],
                mode='lines+markers',
                name='Cluster '+str(j)
            )
            data.append(trace)
        layout = go.Layout(
            title=title,
            xaxis=dict(
                title='Position Along Gene',
                titlefont=dict(
                    size=18,
                    color='#7f7f7f'
                )
            ),
            yaxis=dict(
                title='Probability of Mutation',
                titlefont=dict(
                    size=18,
                    color='#7f7f7f'
                )
            )
        )
        # Plot and generate the html code 
        fig = go.Figure(data=data, layout=layout)
        fig_html = oplot.plot(fig, auto_open=True, output_type='file')
        return fig_html

    def barPlot(self, title="Mutation Distribution Across Positions", ymax=0.5):
        # Compare replicates and generate a plot coloured by bases 
        colors = {'G':'goldenrod', 'U':'mediumturquoise', 'T':'mediumturquoise', 'A':'crimson', 'C':'royalblue', 'N':'black'}
        profile_fig = tools.make_subplots(rows=self.num_clusters, cols=1, shared_yaxes=True)
        log.info(self.num_clusters)
        col_num = 0
        for col in list(self.clusters.columns):
            log.info(col+" data")
            if 'Cluster' in col:
                col_num += 1
                trace1 = go.Bar(
                        x=self.clusters['Position'],
                        y=self.clusters[col],
                        text= self.clusters['Base'],
                        name = col,
                        marker=dict(color=list(self.clusters.Base.map(colors))
                        ),   
                )
                profile_fig.append_trace(trace1,col_num,1)
        profile_fig['layout'].update(title = title)
        for i in range(1,self.num_clusters + 1):
            y_axis = 'yaxis' + str(i)
            profile_fig['layout'][y_axis].update(range=[0,ymax])
        plot_filename = EMClass.split_filename(self.filename).noext + '_plots.html'
        profile_html = oplot.plot(profile_fig , auto_open=False, output_type='file', filename= plot_filename)
        return profile_html

    def mat_barPlot(self, title="Mutation Distribution Across Positions", ymax=0.5):
        # TODO: eps offset 20
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
        colors = {'G':'goldenrod', 'U':'mediumturquoise', 'T':'mediumturquoise', 'A':'crimson', 'C':'royalblue', 'N':'black'}
        # add custom legend for Nucleotide on the first axe
        nuc_handles = [mpatches.Patch(color = color, label = nuc) for nuc, color in colors.items()]
        f, axarr = plt.subplots(self.num_clusters, sharey='row', figsize=(25, 13))
        #f.legend(handles=nuc_handles, loc=1, fontsize=12)
        # if self.num_clusters == 1:
        #     axarr.legend(handles = nuc_handles, loc=1, fontsize=12)
        # else:
        #     axarr[0].legend(handles = nuc_handles, loc=1, fontsize=12)
        log.info(self.num_clusters)
        col_num = 0
        for col in list(self.clusters.columns):
            log.info(col +" data")
            if 'Cluster' in col:
                col_num += 1
                marker = list(self.clusters.Base.map(colors))
                try:
                    cur_ax = axarr[col_num - 1] # current ax
                except TypeError:
                    cur_ax = axarr

                g = sns.barplot(
                    x='Position',
                    y=col,
                    data=self.clusters,
                    palette=marker,
                    ax=cur_ax

                )

                cur_ax.set_ylim(0, ymax)
                # cur_ax.set_xlim(int(list(self.clusters['Position'])[0]), int(list(self.clusters['Position'])[-1]))
                cur_ax.get_xaxis().set_major_locator(ticker.MultipleLocator(10))
                cur_ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter(useOffset=-20))
                plt.setp(cur_ax.get_xaxis().get_offset_text(), visible=False) # Hide the offset text
                g.set_ylabel('Cluster ' + str(col_num), fontsize=16, alpha=0.8)
                sns.set_style('white')  # Set the style to all white background without grids
                g.legend(fontsize = 15)
                # set xlabel only for the bottom sub axe
                if col_num != self.num_clusters:
                    g.set_xlabel('')
                else:
                    g.set_xlabel('Position', fontsize = 15, alpha = 0.8)
                sns.despine()

        plot_filename =  EMClass.split_filename(self.filename).noext + '_plots.eps'
        plt.title(title, y = 3.5, alpha = 0.8)
        plt.savefig(plot_filename, format = 'eps', dpi = 1200)


def run_PrettyPlots(input_clusters_filename, ref=None, sequence=None, ymax=0.5):
    """
    Create a PrettyPlot.
    :param input_clusters_filename:
    :param ref:
    :param sequence:
    :param ymax:
    :return: Void
    """
    clusters = EMClass.clusterGroup(input_clusters_filename)
    cluster_plots = clusterPlots(clusters, ref=ref, sequence=sequence)
    cluster_plots.barPlot(ymax=ymax)
    cluster_plots.mat_barPlot(ymax=ymax)


if __name__ == "__main__":

    import argparse

    '''
    List of arguments and corresponding variables

    Mandatory:
                input_clusters --> input_clusters
                         start --> start
                         end   --> end

    Optional:
                -F, --fasta --> fasta
                  -R, --ref --> ref
    '''

    parser = argparse.ArgumentParser(
        description='Make a pretty plot from a cluster file and a reference.')

    parser.add_argument('input_clusters_filename', type=str,
                        help = 'A cluster file. Example: clusters_Ref_HIV_455_480to900.txt.')

    parser.add_argument('-R','--ref', type=str, dest='ref', default=None,
                        help="Specify if the input_clusters has more than one sequence in its comments.")

    parser.add_argument('-S','--sequence', type=str, dest='sequence', default=None,
                        help="Specify which sequence to use in the plot. "
                             "If None, the same sequence as the clustering will be used.")

    parser.add_argument('-y', '--ymax', type=float, dest='ymax', default=0.5,
                        help="Maximum value for the y-axis.")

    args = parser.parse_args()

    run_PrettyPlots(args.input_clusters_filename, ref=args.ref, sequence=args.sequence, ymax = args.ymax)
