#!/usr/bin/python

import pandas
import plotly
import plotly.graph_objs as go
import numpy as np
import EM_classes as EMClass
import argparse


def scatter_plot(f1, f2=None, cluster_name_1=None, cluster_name_2=None,
                 color_map={('deeppink','prD3A4') : range(709,714)}, output_name=None, seq=None):
    """
    
    :param f1: name of cluster file 1
    :param f2: name of cluster file 2, if None, same as f1
    :param cluster_name_1: name of cluster in f1 header. If None set to 'Cluster_1'
    :param cluster_name_2: name of cluster in f1 header. If None set to 'Cluster_1' if f2!=f1 or 'Cluster_2 if f2==f1
    :param color_map: 
    :param output_name:
    :param seq: Sequence to be used to remove G and Ts. Used instead of sequence found in header.
    :return: 
    """
    if f2 is None:
        f2 = f1
    split_f1 = EMClass.split_filename(f1)
    split_f2 = EMClass.split_filename(f2)
    if cluster_name_1 is None:
        cluster_name_1 = "Cluster_1"
    if cluster_name_2 is None:
        if f2 != f1:
            cluster_name_2 = "Cluster_1"
        else:
            cluster_name_2 = "Cluster_2"
    cluster1 = EMClass.singleCluster(f1, cluster_name_1)
    cluster2 = EMClass.singleCluster(f2, cluster_name_2)

    if seq is None:
        seq = cluster1.get_sequence()
        if not cluster2.get_sequence() == seq:
            print "The 2 clusters have different sequences: \n\n%s \n\n %s" % (seq, cluster2.probs.shape[0])
    if not cluster1.probs.shape[0] == cluster2.probs.shape[0]:
        print "The 2 clusters have different lengths: %s vs %s so the interscetion will be compared" % (cluster1.probs.shape[0], cluster2.probs.shape[0])
    if not (cluster1.probs.index == cluster1.probs.index).all():
        print "Warning: the 2 cluster you are comparing don't have the same indices"
    seq = list(seq.strip())
    # assert(len(seq) == cluster1.probs.shape[0]), "The sequence length (%s) is different from the cluster length (%s)"\
    #                                              % (len(seq), cluster1.probs.shape[0])
    # get the index of intersection
    idx = cluster1.probs.index.intersection(cluster2.probs.index)
    assert list(idx) != [], "There is not overlapping section between the 2 clusters."
    data_all = pandas.concat([cluster1.probs[idx], cluster2.probs[idx]], axis=1)
    data_all.columns = range(data_all.shape[1])
    # truncate the sequence for the intersection
    seq_series = pandas.Series(list(seq))
    seq_series.index = cluster1.probs.index
    cond = [x == 'A' or x == 'a' or x == 'C' or x == 'c' for x in seq_series[idx]]
    data_all = data_all[cond]

    if output_name == None:
        try:
            output_name = 'scatterplot_' + split_f1.path.split('/')[-2] + '_' + cluster1.name + '_VS_' + split_f2.path.split('/')[-2] + '_' + \
                          cluster2.name + '.html'
        except:
            output_name = 'scatterplot_' + split_f1.basename + '_' + cluster1.name + '_VS_' + split_f2.basename + '_' + \
                          cluster2.name + '.html'

    data = []

    ignore_pts = []
    for (color, name) in color_map:
        ignore_pts += color_map[(color, name)]

    inds = []
    ignore_set = set(ignore_pts)
    for val in data_all.index:
        if val not in ignore_set:
            inds.append(val)

    color_map[('grey', 'All Bases')] = inds#list(set(data_all.index) - set(ignore_pts))

    for (color, name) in color_map:
        special_vals = color_map[(color, name)]
        df = data_all.loc[data_all.index.isin(special_vals)]
        x = df[0]
        y = df[1]

        if name == 'All Bases':
            A = np.vstack([x, np.ones(len(x))]).T
            [slope, intercept] = np.linalg.lstsq(A, y)[0]

            best_fit_line = x * slope + intercept
            all_corr = x.corr(y)
            # get R square
            model, resid = np.linalg.lstsq(A, y)[:2]
            r2 = 1 - resid / (y.size * y.var())
            # plot linear regression
            best_fit_trace = go.Scatter(
                x=x,
                y=best_fit_line,
                mode='lines',
                name='Best Fit',
                marker=go.Marker(color='lightgrey'),
            )

            data.append(best_fit_trace)

        trace = go.Scatter(
            x=x,
            y=y,
            text=special_vals,
            mode='markers',
            name=name,
            marker=dict(size=15, color=color, line=dict(width=1,)))
        data.append(trace)

    title = "Regression between Replicates - R^2: %0.4f " % (r2)

    layout = go.Layout(
        title=title,
        xaxis=dict(
            title = split_f1.basename + " " + cluster1.name,
        ),
        yaxis=dict(
            title=split_f2.basename + " " + cluster2.name,
        ),
        showlegend=False
    )

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=output_name, auto_open=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Make a scatter plot out of 2 cluster_probs')

    parser.add_argument('f1',type=str,
        help='String. Cluster file containing the first cluster.')

    parser.add_argument('-f2', '--input_filename_2', type=str, dest='f2', default=None,
                        help='String. Cluster file containing the second cluster. If None, will use f1.')

    parser.add_argument('-c1', '--cluster_name_1', type=str, dest='cluster_name_1', default=None,
                        help='String. Name of the cluster of interest in f1 header. If None, will use "Cluster_1".')

    parser.add_argument('-c2', '--cluster_name_2', type=str, dest='cluster_name_2', default=None,
                        help='String. Name of the cluster of interest in f1 header. If None, will use "Cluster_1" '
                             'unless f2 is the same as f1 in which case it will default to "Cluster_2".')

    parser.add_argument('-S', '--sequence', type=str, dest='seq', default=None,
                        help="String. Sequence to use to remove G and Ts. Specify only if using old cluster files that"
                             "don't have header. This option will overwrite the sequence found in headers.")


    args = parser.parse_args()

    scatter_plot(args.f1, f2=args.f2, cluster_name_1=args.cluster_name_1, cluster_name_2=args.cluster_name_2,
                 seq=args.seq)






