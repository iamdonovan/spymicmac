#!/usr/bin/env python
import argparse
from spymicmac import micmac


def _argparser():
    parser = argparse.ArgumentParser(description="Create MicMac-LocalChantierDescripteur.xml, "
                                                 "AutoCal_Foc-304800_KH9MC.xml, id_fiducial.txt, MeasuresCamera.xml, "
                                                 "and Tmp-SL-Glob.xml files for KH-9 Hexagon Mapping Camera images.",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-joined', action='store_true',
                        help='generate files for joined scene (220x460 mm) instead of half (220x230mm)')
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()
    micmac.init_autocal()
    micmac.create_localchantier_xml()
    micmac.generate_measures_files(joined=args.joined)


if __name__ == "__main__":
    main()
