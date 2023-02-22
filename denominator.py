def smart_denominator(mu, dmuts=3):
    # Computes the actual denominator in polynomial time without generation xA based on recursion
    dens = {}
    res = den_helper(mu, dens, dmuts=dmuts)
    return res, dens


def den_helper(mu, dens, dmuts=3):
    if len(mu) in dens:
        return dens[len(mu)]
    elif len(mu) == 0:
        return 1
    else:
        s1 = (1 - mu[0]) * den_helper(mu[1:], dens, dmuts=dmuts)
        s2 = mu[0] * (1.0 - mu[1:dmuts+1]).prod() * den_helper(mu[dmuts+1:], dens, dmuts=dmuts)
        dens[len(mu)] = s1 + s2
        return dens[len(mu)]


def den_alpha(mu, mu_rev, alpha, dens, dens_rev, dmuts=3):
    return mu[alpha] * den_alpha_side(mu, alpha, dens, dmuts=dmuts) * den_alpha_side(mu_rev, len(mu) - alpha - 1, dens_rev, dmuts=dmuts)


def den_alpha_side(mu, alpha, dens, dmuts=3):
    if len(mu[alpha + dmuts+1:]) == 0:
        return (1.0 - mu[alpha + 1:alpha + dmuts+1]).prod()
    return (1.0 - mu[alpha + 1:alpha + dmuts+1]).prod() * dens[len(mu[alpha + dmuts+1:])]

