CXX ?= g++
FLAGS = -Wall -fopenmp -std=c++20
SRC_DIR = source/
BUILD_DIR = build/
SRC = ${wildcard $(SRC_DIR)*.cpp}
OBJ = ${patsubst %.cpp, $(BUILD_DIR)%.o, $(SRC)}
INCLUDE = -I./include

DRIVERS_DIR = drivers/
DRIVERS_SRC = ${wildcard $(DRIVERS_DIR)*.cpp}
DRIVERS_OUT = ${patsubst %.cpp, %, $(DRIVERS_SRC)}

PAPER_DIR = paper/
PAPER_SRC = ${wildcard $(PAPER_DIR)*.cpp}
PAPER_OUT = ${patsubst %.cpp, %, $(PAPER_SRC)}
PAPER_LIBS = -larmadillo

LIB = fdwh
LIBF = lib$(LIB).a

waveholtz debug: $(OBJ)
	ar rcs $(LIBF) $(OBJ)

debug: FLAGS += -g -DWH_DEBUG

waveholtz: FLAGS += -O3

build/%.o: %.cpp
	@mkdir -p ${dir $@}
	$(CXX) $(FLAGS) -o $@ $< -c $(INCLUDE)

drivers: $(DRIVERS_OUT)
	@mkdir -p solution

paper: $(PAPER_OUT)

$(DRIVERS_DIR)%: $(DRIVERS_DIR)%.cpp $(LIBF)
	$(CXX) $(FLAGS) -o $@ $< $(INCLUDE) -L. -l$(LIB)

$(PAPER_DIR)%: $(PAPER_DIR)%.cpp $(LIBF)
	$(CXX) $(FLAGS) -o $@ $< $(INCLUDE) $(PAPER_LIBS) -L. -l$(LIB)

clean:
	rm -rf build
	rm -f $(LIBF)
	rm -f $(DRIVERS_OUT)