
def calc_percentiles(cnts_dict, percentiles_to_calc=range(101)):
    """Returns [(percentile, value)] with nearest rank percentiles.
    Percentile 0: <min_value>, 100: <max_value>.
    cnts_dict: { <value>: <count> }
    percentiles_to_calc: iterable for percentiles to calculate; 0 <= ~ <= 100
    """
    assert all(0 <= p <= 100 for p in percentiles_to_calc)
    percentiles = []
    num = sum(cnts_dict.values())
    cnts = sorted(cnts_dict.items())
    curr_cnts_pos = 0  # current position in cnts
    curr_pos = cnts[0][1]  # sum of freqs up to current_cnts_pos
    for p in sorted(percentiles_to_calc):
        if p < 100:
            percentile_pos = p / 100.0 * num
            while curr_pos <= percentile_pos and curr_cnts_pos < len(cnts):
                curr_cnts_pos += 1
                curr_pos += cnts[curr_cnts_pos][1]
            percentiles.append((p, cnts[curr_cnts_pos][0]))
        else:
            percentiles.append((p, cnts[-1][0]))  # we could add a small value
    return percentiles
