import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate a water rocket launch.")
    parser.add_argument("word", type=str, help="a word to print")

    args = parser.parse_args()
    print(f"Hello! You said: {args.word}")
