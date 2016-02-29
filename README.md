# TreeFig
Converter that takes a figure of a phylogenetic tree and calculates a tree notation.

## Requirements:
Requires the `xmltodict` and `svg.path` packages.

## Testing:
There are some test .svg images in the test directory. Run the test script to generate their output files.
	```bash test/run_tests.sh```

## Running:
I recommend using ImageMagick's suite of command-line tools to convert any screencap of a tree into a portable bitmap image.

Use the threshold function to render the image down to black and white; I recommend 50%, but you can tweak as needed.
    ```convert input_image -threshold 50% raw.pbm```
    
Then use potrace to vectorize the image:
    ```potrace -s -k 0.8 -W 10 -H 10 -o my_tree.svg raw.pbm```

Run the script with the vectorized file as the argument:
	```python svg_converter.py my_tree.svg```

The output files include `my_tree_raw.svg`, which contains the intermediate steps used in the calculation (this is good for debugging), `my_tree.nex`, and `my_tree.nexml`.

## Caveats:
* As of Feb 28, 2016, the script only processes right-facing rectangular cladograms. 
* It does not OCR the taxon names yet, but it does number them from the top down for easy matching later.
* The script does not handle very short nodes well.
