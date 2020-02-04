import argparse
from shapely.geometry import LineString
from elements import elements_dict
from itertools import combinations
from scipy.integrate import simps
import operator
from pypif import pif
from pypif.obj import *


def find_min_energy_overlap(lines, intersect_points):
    """
    Finds the minimum energy line given the original lines for a defect and the intersection points of those lines

    Args:
        lines (array of arrays): Energy lines for a particular defect. Each line corresponds to a different defect charge.
        intersect_points (array): Points of intersection for overlapping energy lines.

    Returns:
        [min_y_values, min_x_values]: lowest energy line based on area under curve.
    """

    # set x1 and x2, x1 = 0, x2 = bg
    x1 = lines[0][0][0]
    x2 = lines[0][1][0]

    # set to first y value
    y1 = lines[0][0][1]
    y2 = lines[0][1][1]

    # compare y against all other y values
    for pair in lines:
        if pair[0][1] < y1:
            y1 = pair[0][1]
        if pair[1][1] < y2:
            y2 = pair[1][1]

    # create combo of intersect points to generate candidate lines. Assumption is that lowest energy line will have int_points = len(lines)-1
    int_combo = list(combinations(intersect_points, len(lines)-1))

    y_values = [y1]
    x_values = [x1]
    min_y_values = []
    min_x_values = []
    min_area = 100000   # initialize to arb. large number

    # iterate through candidate points to find curve with min area.
    for candidate_points in int_combo:
        for xy in candidate_points:
            y_values.append(xy[1])
            x_values.append(xy[0])
        y_values.append(y2)
        x_values.append(x2)
        area_under_curve = simps(y_values, x_values)

        if area_under_curve < min_area:
            min_area = area_under_curve
            min_y_values = y_values
            min_x_values = x_values

        print("Y: ", y_values, "  X: ", x_values, " | ", round(area_under_curve, 4))

        y_values = [y1]
        x_values = [x1]

    return [min_y_values, min_x_values]


def calculate_intersect_points(lines):
    """
    Calculates intersection points for defects with >1 line/charge state.

    Args:
        lines (array of arrays): Energy lines for a particular defect. Each line corresponds to a different defect charge.

    Returns:
        intersection_points (list): returns a list containing any point where two lines intersect
    """

    intersection_points = []

    for pair in list(combinations(lines, 2)):
        line1 = LineString(pair[0])
        line2 = LineString(pair[1])
        intersect_point = line1.intersection(line2)
        try:
            intersection_points.append([round(intersect_point.x, 3), round(intersect_point.y, 3)])
        except AttributeError:
            print("NO INTERSECTION POINTS")

    return intersection_points


def calc_defect_enthalpy(enthalpies_at_corner, enthalpy_at_mu_0, defect_type, defect_site):
    """
    Calculates defect enthalpy.

        Args:
            enthalpies_at_corner (dict): enthalpies at specific corner given from user/template
            enthalpy_at_mu_0 (str): relative enthalpy at $/mu_0$
            defect_type (str): type of defect. I = interstital, V = vacancy, else = defect is element
            defect_site (str): site of defect. usually element

        Returns:
            defect_enthalpy (float): calc enthalpy of defect
    """

    enthaply_of_site = enthalpies_at_corner[defect_site]

    if defect_type == "I":
        defect_enthalpy = float(enthalpy_at_mu_0)+float(enthaply_of_site)*-1

    elif defect_type == "V":
        defect_enthalpy = float(enthalpy_at_mu_0)+float(enthaply_of_site)

    else:
        defect_enthalpy = float(enthalpy_at_mu_0)+1*float(enthalpies_at_corner[defect_site])+(-1)*float(enthalpies_at_corner[defect_type])

    return round(defect_enthalpy, 4)


def get_values(defect_template):
    """
    Opens template and returns dict of two-column data.
        Args:
            defect_template (closed_file): closed csv file containing user submitted template

        Returns:
            entries (dict): key value pairs for two-column data from template
    """


    entries = {}

    for line in open(defect_template, "r").readlines():
        line = line.strip("\r\n")
        if line.split(",")[0] != "":
            entries[line.split(",")[0]] = line.split(",")[1]

    return entries


def parse_template(defect_template):
    """
    Main parsing function. Produces pifs.

        Args:
            defect_template (closed_file): closed csv file containing user submitted template

        Returns:
            sytems (pifs): list of pifs created from template
    """

    systems = []
    atoms = []
    corners = []
    enthalpies = []
    entries = get_values(defect_template)

    band_gap = float(entries['bg'])

    for k, v in entries.items():

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

    count = 1
    print("NUMBER OF CORNERS: ", len(enthalpies))

    for corner in enthalpies:
        print("\n=====")
        system = ChemicalSystem()
        system.chemical_formula = "".join(atoms)
        system.properties = []
        system.ids = []
        system.ids.append(Id(name="Corner", value=count))
        system.ids.append(Id(name="Corner", value=max(corner.items(), key=operator.itemgetter(1))[0]+"-rich"))
        count += 1
        print("CORNER:", corner, pif.dumps(system.ids))

        # initialize dict. k=defect, v=list of energy values for that defect len == number of charges that defect can take
        unique_defects = {}

        for k, v in entries.items():
            if len(k.split("_")) > 3:
                defect_type = k.split("_")[0]
                site = k.split("_")[1]
                charge = k.split("_")[2]
                index = k.split("_")[3]

                y1_enthalpy_at_0 = float(calc_defect_enthalpy(corner, v, defect_type, site))
                y2_enthalpy_at_ef = round(float(charge)*float(band_gap)+float(y1_enthalpy_at_0), 4)

                try:
                    unique_defects[defect_type+"_"+site].append([[0, y1_enthalpy_at_0], [band_gap, y2_enthalpy_at_ef]])
                except KeyError:
                    unique_defects[defect_type+"_"+site] = [[[0, y1_enthalpy_at_0], [band_gap, y2_enthalpy_at_ef]]]

                print("Defect key: ", k, " Enthalpy at x = 0: ", y1_enthalpy_at_0, " Enthalpy at x = band gap: ", y2_enthalpy_at_ef)

                # create properties and append to system
                system.properties.append(Property(name="$\Delta$H", scalars=[y1_enthalpy_at_0, y2_enthalpy_at_ef], conditions=[Value(name="E$_F$", scalars=[0, band_gap])]))
                defect_enthalpy_prop = Property(name="Defect Enthalpy", scalars=y1_enthalpy_at_0)
                defect_enthalpy_prop.conditions = []
                defect_enthalpy_prop.conditions.append(Value(name="Defect type", scalars=defect_type))
                defect_enthalpy_prop.conditions.append(Value(name="Defect site", scalars=site))
                defect_enthalpy_prop.conditions.append(Value(name="Defect charge", scalars=charge))
                defect_enthalpy_prop.conditions.append(Value(name="Defect index", scalars=index))
                defect_enthalpy_prop.conditions.append(Value(name="Defect label", scalars=k))
                system.properties.append(defect_enthalpy_prop)

        # calc intersection points for overlapping lines
        for k, v in unique_defects.items():
            print("\n-----CALCULATING INTERSECTION POINTS-----")
            print("Defect: ", k, " Number of charge states: ", len(v))
            print("Defect curves: ", v)

            if len(v) >= 2:
                intersection_points = calculate_intersect_points(v)
                print("INTERSECTION POINTS: ", intersection_points)
                low_energy_line = find_min_energy_overlap(v, intersection_points)
                print("LOWEST ENERGY LINE: ", low_energy_line)
                system.properties.append(Property(name="$\Delta$H_2", scalars=low_energy_line[0], conditions=[Value(name="E$_F$_2", scalars=low_energy_line[1])]))
            else:
                print("LOWEST ENERGY LINE: ", [[v[0][0][0], v[0][1][0]], [v[0][0][1], v[0][1][1]]])
                system.properties.append(Property(name="$\Delta$H_2", scalars=[v[0][0][1], v[0][1][1]], conditions=[Value(name="E$_F$_2", scalars=[v[0][0][0], v[0][1][0]])]))

        systems.append(system)
        print("=====")

    return systems


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", nargs="*", help="path to template file")

    args = parser.parse_args()

    for f in args.csv:
        pifs = parse_template(f)
        outfile = f.replace(".csv", ".json")
        pif.dump(pifs, open(outfile, "w"))
        print("PIF DUMPED: ", outfile)

