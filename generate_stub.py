import os
import sys
import argparse
import importlib
from pathlib import Path
from nanobind.stubgen import StubGen

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--module-name", required=True, help="Name of the module to import")
    parser.add_argument("--pyd-path", required=True, help="Path to the directory containing the .pyd/.so file")
    parser.add_argument("--output-file", required=True, help="Path of the top-level .pyi file (others go to same dir)")
    parser.add_argument("--dll-path", help="[Optional] Extra path to search for DLLs (needed on Windows+MinGW)")
    args = parser.parse_args()

    # DLL search path (Windows)
    if args.dll_path:
        if sys.platform == "win32" and sys.version_info >= (3, 8):
            if os.path.isdir(args.dll_path):
                os.add_dll_directory(args.dll_path)

    sys.path.insert(0, args.pyd_path)
    try:
        module = importlib.import_module(args.module_name)
    except ImportError as e:
        print(f"\nFATAL: Failed to import '{args.module_name}'")
        print(f"Full Error: {e}")
        sys.exit(1)

    output_file = Path(args.output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    sg = StubGen(
        module=module,
        recursive=True,
        include_docstrings=True,
        quiet=False,
        output_file=output_file
    )

    sg.put(module)

    with open(output_file, "w", encoding="utf-8") as f:
        f.write(sg.get())

if __name__ == "__main__":
    main()
