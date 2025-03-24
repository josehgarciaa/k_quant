import os

def merge_files(input_directory, output_file):
    """Reads all text-based files from a directory and merges them into a single file."""
    with open(output_file, 'w', encoding='utf-8') as outfile:
        for root, _, files in os.walk(input_directory):
            for file in files:
                file_path = os.path.join(root, file)
                if file.endswith(('.py',)):  
                    try:
                        with open(file_path, 'r', encoding='utf-8') as infile:
                            outfile.write(f"\n\n--- File: {file_path} ---\n\n")
                            outfile.write(infile.read() + "\n")
                    except Exception as e:
                        print(f"Skipping {file_path} due to error: {e}")

if __name__ == "__main__":
    directory = input("Enter the directory path: ")
    output = "merged_output.txt"
    merge_files(directory, output)
    print(f"All files have been merged into {output}")
