CXX ?= g++
FLAGS = -Wall -fopenmp
SRC_DIR = source/
BUILD_DIR = build/
SRC = ${wildcard $(SRC_DIR)*.cpp}
OBJ = ${patsubst %.cpp, $(BUILD_DIR)%.o, $(SRC)}
INCLUDE = -I./include

DRIVERS_DIR = drivers/
DRIVERS_SRC = ${wildcard $(DRIVERS_DIR)*.cpp}
DRIVERS_OUT = ${patsubst %.cpp, %, $(DRIVERS_SRC)}

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

$(DRIVERS_DIR)%: $(DRIVERS_DIR)%.cpp $(LIBF)
	$(CXX) $(FLAGS) -o $@ $< $(INCLUDE) -L. -l$(LIB)

clean:
	rm -rf build
	rm -f $(LIBF)
	rm -f $(DRIVERS_OUT)