EXE_NAME = cpu_fdeb 

CC = mpic++
CXX = $(CC)


CFLAGS += -g
#CFLAGS += -O3
CFLAGS += -Wall
CFLAGS += -DDEBUG
#CFLAGS += -DTIMING

OBJ_DIR = .
CFILES = main.cpp gldraw.cpp point.cpp line.cpp

OBJECT_FILES = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(CFILES))

all: $(OBJECT_FILES)
	$(CXX) -g -o ${EXE_NAME} $(OBJECT_FILES) $(LIB) $(LDFLAGS)

.depend: $(CFLIES)
	$(CXX) -MM -E $(CFLAGS) $(INCLUDE) $(CFILES) > .depend	
	
$(OBJ_DIR)/%.o: %.cpp
	$(CC) -c -o $@ $<

#include .depend

clean:
	rm -f core.* $(OBJ_DIR)/*.o $(EXE_NAME) .depend
