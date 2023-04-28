#!/usr/bin/env python
import argparse
from spymicmac.orientation import block_orientation


def _argparser():
    parser = argparse.ArgumentParser(description="Combine GCPs, Measures files, and Ori directories from multiple "
                                                 "sub-blocks into a single file and orientation.",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('blocks', action='store', type=int, nargs='+', help='the block numbers to combine')
    parser.add_argument('-meas_out', action='store', type=str, default='AutoMeasures',
                        help='the output filename for the Measures file (no extension). (default: AutoMeasures)')
    parser.add_argument('-gcp_out', action='store', type=str, default='AutoGCPs',
                        help='the output filename for the GCP file (no extension). (default: AutoGCPs)')
    parser.add_argument('-fn_mes', action='store', type=str, default='AutoMeasures_block',
                        help='the name pattern of the measures files to combine (default: AutoMeasures_block)')
    parser.add_argument('-fn_gcp', action='store', type=str, default='AutoGCPs_block',
                        help='the name pattern of the GCP files to combine (default: AutoGCPs_block)')
    parser.add_argument('-d', '--dirname', action='store', type=str, default='auto_gcps',
                        help='the output directory where the files are saved (default: auto_gcps)')
    parser.add_argument('-r', '--rel_ori', action='store', type=str, default='Relative',
                        help='the relative orientation to input to GCPBascule (default: Relative -> Ori-Relative)')
    parser.add_argument('-o', '--outori', action='store', type=str, default='TerrainFinal',
                        help='the output orientation from Campari (default: TerrainFinal -> Ori-TerrainFinal)')
    parser.add_argument('-homol', action='store', type=str, default='Homol',
                        help='the Homologue directory to use (default: Homol)')
    parser.add_argument('-ref_dx', action='store', type=float, default=15,
                        help='the pixel resolution of the reference image, in meters. (default: 15)')
    parser.add_argument('-ortho_res', action='store', type=float, default=8,
                        help='the pixel resolution of the orthoimage being used, in meters. (default: 8)')
    parser.add_argument('-allfree', action='store_false',
                        help='run Campari with AllFree set to False')
    parser.add_argument('-max_iter', action='store', type=int, default=1,
                        help='the maximum number of iterations to run. (default: 1)')
    return parser


def main():
    parser = _argparser()
    args = parser.parse_args()

    block_orientation(**vars(args))


if __name__ == "__main__":
    main()
