def find_p_value(number_EQ,probability,value):
    """
        Calculate the p value (one-sided) base on Bernoulli
    """
    from math import erf,sqrt

    if value>number_EQ*probability:
        p_value = 0.5-0.5*erf((value-number_EQ*probability)/sqrt(2*number_EQ*probability*(1-probability)))
    else:
        p_value = 0.5+0.5*erf((value-number_EQ*probability)/sqrt(2*number_EQ*probability*(1-probability)))




    return p_value
