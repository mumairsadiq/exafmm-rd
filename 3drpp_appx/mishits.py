import sys

def find_missing_lines(file1_path, file2_path):
    try:
        with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
            file1_lines = set(file1.readlines())
            file2_lines = set(file2.readlines())

        missing_lines = file1_lines - file2_lines

        print(f"Lines in {file1_path} but not in {file2_path}:")
        for line in sorted(missing_lines):  # Sorting to keep the output consistent
            print(line.strip())

    except FileNotFoundError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <file1_path> <file2_path>")
        sys.exit(1)

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]

    find_missing_lines(file1_path, file2_path)
