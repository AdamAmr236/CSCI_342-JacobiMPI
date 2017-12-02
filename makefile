# Environment
# "development": build with debugging symbols and minimal optimization
# "production": build without debugging symbols and maximum optimization
#environment ?= development

# Source files directory
SRCDIR := src
# Header files directory
INCDIR := inc
# Object files directory
OBJDIR := obj

# Executable
EXE := jacobi
# Source files
SRC := $(notdir $(shell find $(SRCDIR) -name '*.c'))
# Object files
OBJ := $(addprefix $(OBJDIR)/, $(SRC:.c=.o))
# Dependeny files
DEP := $(OBJ:.o=.d)

# Prerequisites path
vpath %.c $(shell find $(SRCDIR) -type d)
vpath %.h $(INCDIR)

# Compiler
CC = gcc
# Compiler flags
CFLAGS := -fopenmp -c -pedantic -std=gnu11 -m64 $(addprefix -I, $(INCDIR))
# Environment specific compiler flags
ifeq ($(environment), development)
	CFLAGS += -Wall
	CFLAGS += -Wextra
	CFLAGS += -Winline
	CFLAGS += -Wmissing-prototypes
	CFLAGS += -Wnested-externs
	CFLAGS += -Wredundant-decls
	CFLAGS += -Wshadow
	CFLAGS += -Wstrict-prototypes
	CFLAGS += -Wundef
	CFLAGS += -Wunreachable-code
	CFLAGS += -ggdb3
	CFLAGS += -D_DEBUG
else
	CFLAGS += -flto
	CFLAGS += -Ofast
	CFLAGS += -march=core-avx-i
	CFLAGS += -fopt-info-vec=optimizations.opt
endif
# Output option flags
OUTPUT_OPTION = -MMD -MP -o $@

LD = $(CC)
# Linker flags
LDFLAGS := -m64 -fopenmp 
# Environment specific linker flags
ifeq ($(environment), development)
	LDFLAGS += -ggdb3
else
	LDFLAGS += -flto
	LDFLAGS += -Ofast
	LDFLAGS += -march=core-avx-i
endif
# Libraries
LIBS = m
# Linker libraries
LDLIBS := $(addprefix -l, $(LIBS))

# Compile command
COMPILE = $(CC) $(CFLAGS) $(OUTPUT_OPTION)
# Link command
LINK = $(LD) -o $@

# Non-generating rules
.PHONY : all clean cleanexe cleanobj

# Default rule
all : $(EXE)

# Link and build executable
$(EXE) : $(OBJ)
	$(LINK) $^ $(LDFLAGS) $(LDLIBS)

# Create objects directory
$(OBJDIR):
	mkdir -p $@

# Bootstrap dependency generation
$(OBJDIR)/%.o : %.c | $(OBJDIR)
	$(COMPILE) $<

# Include auto-generated dependencies
# The include must be last to override the previous rules
-include $(DEP)

# Remove all artifacts
clean : cleanexe cleanobj

# Remove executable
cleanexe :
	$(RM) $(EXE)

# Remove objects and dependencies
cleanobj :
	$(RM) $(OBJ) $(DEP)
	$(RM) -r $(OBJDIR)
