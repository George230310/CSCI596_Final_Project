CC = g++
CFLAGS = -c -I./glm -I./png++ `libpng-config --cflags`
LDFLAGS = -lGL -lGLU -lglut -I./glm -I./png++ `libpng-config --ldflags`

all: raytracer

raytracer: raytracer.o Sphere.o
	$(CC) -o raytracer raytracer.o Sphere.o $(LDFLAGS)

raytracer.o: raytracer.cpp raytracer.h Renderable.h
	$(CC) $(CFLAGS) raytracer.cpp

Sphere.o: Sphere.cpp Renderable.h
	$(CC) $(CFLAGS) Sphere.cpp

clean:
	rm -f raytracer raytracer.o Sphere.o
