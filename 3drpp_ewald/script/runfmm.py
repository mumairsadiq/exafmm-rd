import os
import sys
import subprocess
import argparse
import datetime
from rich import print
from rich.markup import escape

def run(cmd : str):
    subprocess.run(
        cmd, 
        shell=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=f'{sys.argv[0]}')
    parser.add_argument('exec', type=str, help='name of execution')
    parser.add_argument('-n', '--num', type=int, default=1000, help='number of particles')
    parser.add_argument('-a', '--algorithm', type=str, default="fed", help='algorithms')
    parser.add_argument('-P', type=int, default=6, help='P of FMM')
    parser.add_argument('-i', '--image', type=int, default=0, help='number of images')
    parser.add_argument('-r', '--rega', type=float, default=0.0, help='size of regularization')
    parser.add_argument('--seed', type=int, help='random seed')
    parser.add_argument('-v', '--verbose', type=int, default=0, help='verbose')
    parser.add_argument('--reg_image0_type', help='type of reg_image0_type')
    args = parser.parse_args()
    
    dt_now = datetime.datetime.now()
    result_path = os.path.abspath("./result")
    result_filename = f"{result_path}/{args.exec}_{dt_now.strftime('%Y-%m-%d-%H-%M-%S')}.log"

    cmd  = f"./build/test/{args.exec}"
    cmd += f" -n {args.num}"
    cmd += f" -a {args.algorithm}"
    cmd += f" -P {args.P}"
    cmd += f" -i {args.image}"
    cmd += f" -r {args.rega}"
    cmd += f" -v {args.verbose}"
    cmd += f" --res_filepath={result_filename}"
    if args.seed:
        cmd += f" --seed {args.seed}"
    if args.reg_image0_type:
        cmd += f" --reg_image0_type={args.reg_image0_type}"
    
    print(f"{cmd = }")
    
    os.makedirs("result", exist_ok=True)
    run(cmd)

    
 