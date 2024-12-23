.SILENT:

CXX = g++
CXXFLAGS = -std=c++11 -I./src/packages -I./src/packages/eigen-3.4.0
LATEX = xelatex
SRC_DIR = src
DOC_DIR = doc
PCH_HPP = $(SRC_DIR)/packages/pch.hpp

all: pch run report
	echo "All Done."

pch: $(SRC_DIR)/packages/pch.hpp
    $(CXX) $(CXXFLAGS) -c $(PCH_HPP) -o $(SRC_DIR)/packages/pch.hpp.gch

run: $(SRC_DIR)/problemA/A \
     $(SRC_DIR)/problemC/C \
     $(SRC_DIR)/problemD/D \
     $(SRC_DIR)/problemE/E \
     $(SRC_DIR)/problemF/F \
     $(SRC_DIR)/requirement5/require5 \
     $(SRC_DIR)/test/test \
     $(SRC_DIR)/intersect/intersect_test
    echo "Running problem A..."
    cd $(SRC_DIR)/problemA && ./A
    cd $(SRC_DIR)/problemA && python plotA.py
    cd $(SRC_DIR)/problemA && python plotA_leastSquare.py
    echo "Running problem C..."
    cd $(SRC_DIR)/problemC && ./C
    cd $(SRC_DIR)/problemC && python plotC.py
    echo "Running problem D..."
    cd $(SRC_DIR)/problemD && ./D
    echo "Running problem E..."
    cd $(SRC_DIR)/problemE && ./E
    cd $(SRC_DIR)/problemE && python plotE.py
    cd $(SRC_DIR)/problemE && python plotE_error.py
    echo "Running problem F..."
    cd $(SRC_DIR)/problemF && ./F
    cd $(SRC_DIR)/problemF && python plotF_trunc.py
    cd $(SRC_DIR)/problemF && python plotF_divide.py
    echo "Running requirement5..."
    cd $(SRC_DIR)/requirement5 && ./require5
    cd $(SRC_DIR)/requirement5 && python plotRequire5.py
    echo "Running test..."
    cd $(SRC_DIR)/test && ./test
    cd $(SRC_DIR)/test && python plotTest.py
    echo "Running intersect test..."
    cd $(SRC_DIR)/intersect && ./intersect_test
    cd $(SRC_DIR)/intersect && python plotIntersect.py
	echo "Programm Done."


$(SRC_DIR)/problemA/A: $(SRC_DIR)/problemA/A.cpp $(SRC_DIR)/packages/pch.hpp.gch
    echo "Compiling problem A..."
    $(CXX) $(CXXFLAGS) -include $(PCH_HPP) $< -o $@

$(SRC_DIR)/problemC/C: $(SRC_DIR)/problemC/C.cpp $(SRC_DIR)/packages/pch.hpp.gch
    echo "Compiling problem C..."
    $(CXX) $(CXXFLAGS) -include $(PCH_HPP) $< -o $@

$(SRC_DIR)/problemD/D: $(SRC_DIR)/problemD/D.cpp $(SRC_DIR)/packages/pch.hpp.gch
    echo "Compiling problem D..."
    $(CXX) $(CXXFLAGS) -include $(PCH_HPP) $< -o $@

$(SRC_DIR)/problemE/E: $(SRC_DIR)/problemE/E.cpp $(SRC_DIR)/packages/pch.hpp.gch
    echo "Compiling problem E..."
    $(CXX) $(CXXFLAGS) -include $(PCH_HPP) $< -o $@

$(SRC_DIR)/problemF/F: $(SRC_DIR)/problemF/F.cpp $(SRC_DIR)/packages/pch.hpp.gch
    echo "Compiling problem F..."
    $(CXX) $(CXXFLAGS) -include $(PCH_HPP) $< -o $@

$(SRC_DIR)/requirement5/require5: $(SRC_DIR)/requirement5/require5.cpp $(SRC_DIR)/packages/pch.hpp.gch
    echo "Compiling requirement5..."
    $(CXX) $(CXXFLAGS) -include $(PCH_HPP) $< -o $@

$(SRC_DIR)/test/test: $(SRC_DIR)/test/test.cpp $(SRC_DIR)/packages/pch.hpp.gch
    echo "Compiling test..."
    $(CXX) $(CXXFLAGS) -include $(PCH_HPP) $< -o $@

$(SRC_DIR)/intersect/intersect_test: $(SRC_DIR)/intersect/intersect_test.cpp $(SRC_DIR)/packages/pch.hpp.gch
    echo "Compiling intersect test..."
    $(CXX) $(CXXFLAGS) -include $(PCH_HPP) $< -o $@

report:
    echo "Compiling LaTeX report..."
    cd $(DOC_DIR) && $(LATEX) report.tex && $(LATEX) report.tex
    echo "Compiling Design doc..."
    cd $(DOC_DIR) && $(LATEX) design.tex && $(LATEX) design.tex
	echo "Report Done."

clean:
    rm -f figure/problemA/*.png figure/problemC/*.png figure/problemE/*.png figure/requirement5/*.png figure/problemF/*.png figure/test/*.png figure/intersect/*.png
    rm -f output/problemA/*.csv output/problemC/*.csv output/problemD/*.csv output/problemE/*.csv output/requirement5/*.csv output/problemF/*.csv output/test/*.csv output/intersect/*.txt
    rm -f doc/*.aux doc/*.log doc/*.synctex.gz doc/*.pdf
    rm -f src/problemA/A src/problemC/C src/problemD/D src/problemE/E src/problemF/F src/requirement5/require5 src/test/test src/intersect/intersect_test
    rm -f src/packages/pch.hpp.gch
	echo "Cleaned."