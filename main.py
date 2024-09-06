from qsieve import quadratic_sieve
import time
import argparse


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--number", type=int, required=True, help="An integer number to run Quadratic Sieve")
    args = parser.parse_args()

    n = args.number

    start = time.time()
    f1, f2 = quadratic_sieve(n)
    end = time.time()

    print("Time needed in seconds: ", end - start, " to process a number with ", len(str(n)), " digits")

    if f1 == f2 is None:
        print("No solution found for ", n)
    else:
        print("The factors of ", n, "are:", f1, "*", f2)

