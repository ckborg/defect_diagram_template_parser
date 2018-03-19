import argparse
from elements import elements_dict
from pypif import pif
from pypif.obj import *
from citrination_client import PifSystemReturningQuery

def calc_defect_enthalpy(enthalpies_at_corner, enthaply_at_mu_0, defect_type, defect_site):

    enthaply_of_site = enthalpies_at_corner[defect_site]

    if defect_type == "I":
        defect_enthalpy = float(enthaply_at_mu_0)+float(enthaply_of_site)*-1

    elif defect_type == "V":
        defect_enthalpy = float(enthaply_at_mu_0)+float(enthaply_of_site)

    else:
        defect_enthalpy = float(enthaply_at_mu_0)+1*float(enthalpies_at_corner[defect_site])+(-1)*float(enthalpies_at_corner[defect_type])

    return defect_enthalpy


def get_values(defect_template):

    entries = {}

    for line in open(defect_template, "rb").readlines():
        line = line.strip("\r\n")
        if line.split(",")[0] != "":
            entries[line.split(",")[0]] = line.split(",")[1]

    return entries


def parse_template(defect_template):

    systems = []
    atoms = []
    corners = []
    enthalpies = []
    entries = get_values(defect_template)

    band_gap = entries['bg']

    for k, v in entries.iteritems():

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

    count = 0
    for corner in enthalpies:
        print corner
        count += 1
        system = ChemicalSystem()
        system.chemical_formula = "".join(atoms)
        system.properties = []
        system.ids = [Id(name="Corner", value=count)]
        print pif.dumps(system.ids)
        print system.chemical_formula

        print "Values at Corner: ", corner

        for k, v in entries.iteritems():
            if len(k.split("_")) > 3:

                defect_type = k.split("_")[0]
                site = k.split("_")[1]
                charge = k.split("_")[2]
                index = k.split("_")[3]
                y1_enthalpy_at_0 = calc_defect_enthalpy(corner, v, defect_type, site)
                y2_enthalpy_at_ef = float(charge)*float(band_gap)+float(y1_enthalpy_at_0)
                print "DEFECT ENTHALPY: ", k, y1_enthalpy_at_0, y2_enthalpy_at_ef
                system.properties.append(Property(name="Defect enthalpy", scalars=y1_enthalpy_at_0))
                system.properties.append(Property(name="Defect type", scalars=defect_type))
                system.properties.append(Property(name="Defect site", scalars=site))
                system.properties.append(Property(name="Defect charge", scalars=charge))
                system.properties.append(Property(name="Defect index", scalars=index))
                system.properties.append(Property(name="$\delta$H", scalars=[y1_enthalpy_at_0, y2_enthalpy_at_ef], conditions=[Value(name="E$_F$", scalars=[0, band_gap])]))

        systems.append(system)


    return systems


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", nargs="*", help="path to template file")

    args = parser.parse_args()

    for f in args.csv:
        pifs = parse_template(f)
        outfile = f.replace(".csv", ".json")
        pif.dump(pifs, open(outfile, "w"))
        print "PIF DUMPED: ", outfile

