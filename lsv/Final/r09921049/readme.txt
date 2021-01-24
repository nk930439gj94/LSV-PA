Files:
    report.pdf
    download.sh     // command for downloading the source code
    slide.pptx
    papers/


Download:
    ./download.sh


Compile:
    make


Run:
    The implementation is integrated into ABC.
    Command "lsv_esop -h" summerize the usage.
    usage: lsv_esop [-O num] [-h]
    	         esop sythesis
    	-O num : (optional) the 0-based number of the PO to convert to esop [default = the first PO]
    	-p     : toggles print esop [default = no]
    	-h     : print the command usage

    Example:
    read esop_benchmarks/c432.blif
    strash
    dc2
    lsv_esop -O 0 -p

Referenced materials:
    Benchmarks are in "esop_benchmarks" folder.
    Papers:
        Scaling-up ESOP Synthesis for Quantum Compilation
        Logic Synthesis for Quantum Computing
        Pseudo-Kronecker Expressions for Symmetric Functions
        Fast Heuristic Minimization of Exclusive-Sums-of-Products

