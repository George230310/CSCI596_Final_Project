CC = g++
CFLAGS = -c -I./glm -I./png++ `libpng-config --cflags`
LDFLAGS = -lGL -lGLU -lglut -I./glm -I./png++ `libpng-config --ldflags`

all: raytracer

raytracer: raytracer.o
	$(CC) -o raytracer raytracer.o $(LDFLAGS)

raytracer.o: raytracer.cpp
	$(CC) $(CFLAGS) raytracer.cpp

clean:
	rm -f raytracer raytracer.o
