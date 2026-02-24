import os

# Get the directory where the script is located, then go to repo root
# This assumes your script is in /scripts/map_gen.py
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTPUT_FILE = os.path.join(BASE_DIR, "auto_repo_map.md")

IGNORE = {'.git', 'node_modules', '__pycache__', 'venv', 'dist', 'build'}

def generate_map():
    with open(OUTPUT_FILE, "w") as f:
        f.write("# Repository Map\n")
        f.write("> [!NOTE]\n> This file is auto-generated on every commit. Do not edit manually.\n\n")
        
        for root, dirs, files in os.walk(BASE_DIR):
            dirs[:] = [d for d in dirs if d not in IGNORE]
            
            # Calculate depth for indentation
            rel_path = os.path.relpath(root, BASE_DIR)
            if rel_path == ".":
                level = 0
            else:
                level = rel_path.count(os.sep) + 1
            
            indent = '  ' * level
            f.write(f"{indent}- **{os.path.basename(root)}/**\n")
            
            for file in files:
                if file.endswith(('.py', '.js', '.ts', '.go', '.java', '.rb', '.md', ".toml")):
                    f.write(f"{indent}  - {file}\n")
                    # Try to extract signatures
                    try:
                        with open(os.path.join(root, file), 'r', errors='ignore') as code:
                            for line in code:
                                clean = line.strip()
                                # Common patterns for functions/classes
                                if clean.startswith(('def ', 'class ', 'export function', 'interface ', 'async def')):
                                    f.write(f"{indent}    - `{clean.split('{')[0].strip()}`\n")
                    except:
                        pass

if __name__ == "__main__":
    generate_map()