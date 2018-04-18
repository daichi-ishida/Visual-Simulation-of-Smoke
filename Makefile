CXX := g++
BUILD_TYPE := Release

CXX_DEBUG_FLAGS := -g -O0 -Wall
CXX_RELEASE_FLAGS := -s -O2

CXXFLAGS := -MMD -MP -std=gnu++14
BINDIR := bin
SRCDIR := src
OBJDIR := obj
INCLUDE := -I./include

DEBUG := gdb

LIBS := 
SRCS := $(wildcard $(SRCDIR)/*.cpp)
OBJS := $(addprefix $(OBJDIR)/, $(notdir $(SRCS:.cpp=.o)))
DEPS := $(OBJS:.o=.d)

ifeq ($(OS),Windows_NT)
FFTW_PATH := C:/Libraries/fftw-3.3.5-dll64
GLFW_PATH := C:/Libraries/glfw-3.2.1
GLM_PATH := C:/Libraries/glm-0.9.8.5
GLEW_PATH := C:/Libraries/glew-2.1.0
INCLUDE	+= -IC:/MinGW/include -I$(FFTW_PATH) -I$(GLFW_PATH) -I$(GLM_PATH) -I$(GLEW_PATH)
LIBS := -LC:/MinGW/lib -L$(FFTW_PATH) -L$(GLFW_PATH) -L$(GLM_PATH) -L$(GLEW_PATH) -lfftw3-3 -lfftw3f-3 -lfftw3l-3 -lopengl32 -lglew32 -lglfw3
EXECUTABLE	:= main.exe
RM := cmd //C del
else
INCLUDE	+= -I/usr/include/fftw3
LIBS := -lfftw3 -lfftw3f -lfftw3l -lGL -lGLEW -lglfw3 -lX11 -lXxf86vm -lXrandr -lpthread -lXi -ldl -lXinerama -lXcursor
EXECUTABLE	:= main
RM := rm -f
endif

ifeq ($(BUILD_TYPE),Release)
  CXXFLAGS += $(CXX_RELEASE_FLAGS)
else ifeq ($(buildtype),Debug)
  CXXFLAGS += $(CXX_DEBUG_FLAGS)
else
  $(error buildtype must be release, debug, profile or coverage)
endif

.PHONY : all
all: $(BINDIR)/$(EXECUTABLE)

$(BINDIR)/$(EXECUTABLE): $(OBJS) $(LIBS)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

.PHONY : debug
debug :
	$(DEBUG) $(BINDIR)/$(EXECUTABLE)

.PHONY : run
run:
	$(BINDIR)/$(EXECUTABLE)

.PHONY : clean
clean:
	$(RM) $(BINDIR)/$(EXECUTABLE) $(OBJS) $(DEPS)

-include $(DEPS)
