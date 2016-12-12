sources=Utils.cpp Options.cpp LineReader.cpp SNPDistributor.cpp Engine.cpp Search.cpp SearchMI.cpp main.cpp

src_dir = .
objs_dir = objs

#EpistSearch
objs += $(patsubst %.cpp, $(objs_dir)/%.cpp.o, $(sources))

mylibs = -lpthread -lm -lz -lgomp
EXEC = 3SNPInt-Exh

#compile optionsK
CCFLAGS += -O3 -fopenmp
CC = g++

all: dir $(objs) $(objs_wrapper)
	$(CC) $(CCFLAGS) -o $(EXEC) $(objs) $(mylibs)
	strip $(EXEC)

dir:
	mkdir -p $(objs_dir)


clean:
	-rm -rf $(objs_dir) $(EXEC)

$(objs_dir)/%.cpp.o: $(src_dir)/%.cpp
	$(CC) $(CCFLAGS) -o $@ -c $<


