#!/usr/bin/env python
import argparse
import spymicmac.micmac as mt


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
    mt.init_autocal()
    mt.create_localchantier_xml()
    mt.generate_measures_files(joined=args.joined)


if __name__ == "__main__":
    main()
