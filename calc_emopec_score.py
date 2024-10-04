import emopec
import sys
import argparse

def main(args: dict):

    sd = args.sd
    dist = args.dist

    if len(sd) != 6:
        print(f"The length must be SD sequence ({sd}) must be 6", file=sys.stderr)
        exit(0)
        
    print(emopec.get_expression(sd, dist))
    

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('sd', help='sd sequence (6mer)')
    parser.add_argument('dist', type=int, help='distance from ATG')
    args = parser.parse_args()

    main(args)

