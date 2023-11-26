CC = g++
CFLAGS = -Ofast -c -I./glm  -I./libs/libpng/include -I./libs/libpng/lib  -I./png++ 
LDFLAGS = -static-libgcc -static-libstdc++ -I./libs/libpng/include -L./libs/libpng/lib -I./glm -I./png++ -lpng

all: raytracer

raytracer: raytracer.o Sphere.o Triangle.o
	$(CC) -o raytracer raytracer.o Sphere.o Triangle.o $(LDFLAGS)

raytracer.o: raytracer.cpp raytracer.h Renderable.h
	$(CC) $(CFLAGS) raytracer.cpp

Sphere.o: Sphere.cpp Renderable.h
	$(CC) $(CFLAGS) Sphere.cpp

Triangle.o: Triangle.cpp Renderable.h
	$(CC) $(CFLAGS) Triangle.cpp

clean:
	rm -f raytracer raytracer.o Sphere.o Triangle.o
