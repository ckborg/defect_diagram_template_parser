import sys
import pandas as pd
import numpy as np
from pypif import pif
from pypif.obj import *

def calculate_dfe(ps, header_dict, data_dfe, eg):
    """function that return defect formation energy as function of fermi-energy for each defect at given phase stability point

    input: ps (integer) "phase stability point"
    returns: x (float) "fermi-energy"
             y (float) "defect formation energy"
             name (string) "defect label"
    """

    # rearranging data
    name = []
    d = np.zeros([len(data_dfe), 3])

    name.append(data_dfe.values[0][0])

    d[0, 0] = 0
    d[0, 1] = data_dfe.values[0][2]
    d[0, 2] = data_dfe.values[0][3] + data_dfe.values[0][4] * header_dict['dmuA'][ps] + data_dfe.values[0][5] * \
              header_dict['dmuB'][ps] + data_dfe.values[0][6] * header_dict['dmuC'][ps] + data_dfe.values[0][7] * \
              header_dict['dmuD'][ps] + data_dfe.values[0][8] * header_dict['dmuE'][ps] + data_dfe.values[0][9] * \
              header_dict['dmuF'][ps] + data_dfe.values[0][10] * header_dict['dmuG'][ps]

    count = 0

    for i in range(1, len(data_dfe)):
        if data_dfe.values[i][0] != data_dfe.values[i - 1][0] or data_dfe.values[i][1] != data_dfe.values[i - 1][1]:
            count = count + 1
            name.append(data_dfe.values[i][0])

        d[i, 0] = count
        d[i, 1] = data_dfe.values[i][2]
        d[i, 2] = data_dfe.values[i][3] + data_dfe.values[i][4] * header_dict['dmuA'][ps] + data_dfe.values[i][5] * \
                  header_dict['dmuB'][ps] + data_dfe.values[i][6] * header_dict['dmuC'][ps] + data_dfe.values[i][7] * \
                  header_dict['dmuD'][ps] + data_dfe.values[i][8] * header_dict['dmuE'][ps] + data_dfe.values[i][9] * \
                  header_dict['dmuF'][ps] + data_dfe.values[i][10] * header_dict['dmuG'][ps]

    unq = np.unique(d[:, 0])

    # storing dH vs. EF for each unique defect
    x = np.linspace(0, eg, 100)
    y = np.zeros([len(x), len(unq)])

    for kk in range(len(unq)):
        dummy = np.full((len(x)), 100.)
        for jj in range(len(d)):
            if d[jj, 0] == unq[kk]:
                for ii in range(len(x)):
                    val = d[jj][1] * x[ii] + d[jj][2]
                    if val <= dummy[ii]: dummy[ii] = val

        y[:, kk] = dummy[:]

    # function returns x, y, name
    return (x, y, name)


def parse_template(f):
    # read csv
    print(f)
    data = pd.read_csv(f)
    df1 = pd.DataFrame(data)

    # find last row of DFE data
    for i in range(len(data.defect)):
        if pd.isnull(data.defect.values[i]) == False:
            end = i

    # find last row of phasestability data
    for j in range(len(data.point)):
        if pd.isnull(data.point.values[j]) == False:
            endp = j

    # read dfe at dmu[i]=0
    data_dfe = pd.concat(
        [df1['defect'][:end + 1], df1['site'][:end + 1], df1['charge'][:end + 1], df1['dHvbm'][:end + 1],
         df1['na'][:end + 1], df1['nb'][:end + 1], df1['nc'][:end + 1], df1['nd'][:end + 1], df1['ne'][:end + 1],
         df1['nf'][:end + 1], df1['ng'][:end + 1]], axis=1)

    # bandgap: contained in a metadata file
    eg = data.bandgap.values[0]

    # compound: contained in a metadata file
    compound = data.formula[0]

    # read phase stability information
    data_ps = pd.concat(
        [df1['point'][:endp + 1], df1['dmuA'][:endp + 1], df1['dmuB'][:endp + 1], df1['dmuC'][:endp + 1],
         df1['dmuD'][:endp + 1], df1['dmuE'][:endp + 1], df1['dmuF'][:endp + 1], df1['dmuG'][:endp + 1]], axis=1)

    headers_ps = list(data_ps)

    header_dict = {}
    for h in headers_ps:
        header_dict[h] = data_ps['%s' % h].tolist()


    # choose number of phase stability points to plot
    ps = int(data.phase_stabilty_points.values[0])

    system = ChemicalSystem()
    system.properties = []
    system.chemical_formula = compound

    for stability_point in range(ps):
        x, y, names = calculate_dfe(stability_point, header_dict, data_dfe, eg)

        for i in range(len(y[0, :])):
            prop = Property(name=('$\Delta H_{\mathrm{D,q}}$ at phase stability point '+ str(stability_point)),
                            scalars=[round(float(j), 5) for j in y[:, i]], units='eV')
            cond1 = Value(name='$E_{F}$', scalars=[round(float(k), 5) for k in x], units='eV')
            cond2 = Value(name='Defect label', scalars=names[i])
            prop.conditions = [cond1, cond2]
            system.properties.append(prop)

    return system

system = parse_template(sys.argv[1])
pif.dump(system, open(sys.argv[1].replace('.csv', '.json'), 'w'))