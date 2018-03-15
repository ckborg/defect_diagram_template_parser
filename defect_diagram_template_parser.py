import argparse
from elements import elements_dict


def get_values(defect_template):

    entries = {}

    for line in open(defect_template, "rb").readlines():
        line = line.strip("\r\n")
        if line.split(",")[0] != "":
            entries[line.split(",")[0]] = line.split(",")[1]

    return entries


def parse_template(defect_template):

    vacancy_charge = 1
    int_charge = -1

    atoms = []
    corners = []
    enthalpies = []
    entries = get_values(defect_template)

    band_gap = entries['bg']

    for k,v in entries.iteritems():

        if len(k.split("-")) > 1:
            if k.split("-")[0] in elements_dict.values() and k.split("-")[0] not in atoms:
                atoms.append(k.split("-")[0])
            if k.split("-")[1].isdigit() and k.split("-")[1] not in corners:
                corners.append(k.split("-")[1])

    for corner in corners:
        enthaplies_at_corner = {}
        for atom in atoms:
            enthaplies_at_corner[atom] = entries[atom+"-"+corner]
        enthalpies.append(enthaplies_at_corner)

    for corner in corners:
        print "Values at Corner: ", corner

        for k,v in entries.iteritems():
            if len(k.split("_")) > 3:
                defect_type = k.split("_")[0]
                site = k.split("_")[1]
                charge = k.split("_")[2]
                index = k.split("_")[3]
                enthalpy = v

            print defect_type, site, charge, index, enthalpy




    return []


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", nargs="*", help="path to template file")

    args = parser.parse_args()

    for f in args.csv:
        pifs = parse_template(f)
