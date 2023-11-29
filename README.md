# CSCI596 Final Project: Parallel Ray Tracer

## Authors: Suvi Marathe, Leyu Xu

### Intro: What is Ray Tracing?
Ray tracing is one of the techniques for global illumination in the field of Computer Graphics. In Computer Graphics, global illumination refers to a set of 3D lighting algorithms 
that are capable of not only simulating the light coming directly from a light source but also those reflected off other surfaces in the scene. The concept of global illumination 
exists in contrast to local illumination, in which only the light coming directly from a light source is simulated. In a nutshell, ray tracing algorithm, in its simplest form, allows 
a 3D object in the scene to cast shadow on another 3D object.  
![Global Illumination](readme_images/spheres_casting_soft_shadows_on_each_other_SSAA.png)  
*A 3D Scene Produced by Our Ray Tracer*

### Basic Science of Ray Tracing
The basic science of ray tracing without any mathematics is explained as follows. First of all, we have a very simple 3D scene where we have different objects such as spheres and triangles 
in addition to a point light source. Then, a camera is positioned in the scene where we want our eyes to be looking. Subsequently, let's imagine a plane (or a computer monitor) in front 
of the camera on which the image of the scene will be generated. This plane will be divided to many small grids (or pixels), and through each grid, we will fire a ray from the camera position. 
If this ray does not hit any object in the scene, the grid (pixel) it came from will be colored whatever default color that we assigned. However, if this ray does hit an object, another ray will 
be fired from the hit point towards the point light source. If this second ray reaches the light source without hitting any other objects, then we know that the hit point is not blocked by any 
other objects and its color should be evaluated based on its material property and the color of the light source. Yet, if the second ray hits some other object along the way, then we know that 
the hit point is blocked from the light source, and is therefore not lighted (its corresponding pixel will be completely black with only a single point light source).  
![Ray Tracing Illustration](https://d29g4g2dyqv443.cloudfront.net/sites/default/files/pictures/2018/RayTracing/ray-tracing-image-1.jpg)  
*A Diagram for Ray Tracing Cited from NVIDIA's official website*  

### Our Core Implementations
#### Ray-Sphere Intersection & Ray-Triangle Intersection
At the heart of our ray tracer lies the intersection algorithm between the ray itself and other 3D objects. Since spheres and triangles are some of the most common 3D objects, we are interested in 
developing algorithms that can check ray-sphere and ray-triangle intersections. For the ray-sphere intersection algorithm, the basic idea is very simple. Imagine the center of the sphere being 
a point, and what we want to know is whether the closest distance between any point on the ray and the center of the sphere is less or equal to the radius of the sphere. If it is, then we know that we 
have at least one intersection (usually we have two, one coming into the sphere and one going out). Then, with some more mathematical efforts, we're able to get the exact positions of the intersections 
and our algorithm will choose the closest one to return as the result.  
![Ray-Sphere Intersection](https://www.scratchapixel.com/images/ray-simple-shapes/rayspherecases.png?)  
*A Diagram for Ray Sphere Intersection Cited from Scratchapixel 3.0*  

When it comes to the ray-triangle intersection algorithm, things would get a bit more complicated. First of all, the algorithm will calculate if the ray intersect with the plane on which the triangle 
resides. If it does, then the algorithm will calculate the **barycentric coordinates** of intersection point. A set of barycentric coordinates expresses the intersection point as a combination of the three 
vertices of the triangle using the subarea ratios of the sub-triangles formed by the intersection point and any two vertices.  
![Ray-Triangle Intersection](https://la.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/49670/versions/3/screenshot.jpg)  
*A Diagram for Ray Triangle Intersection Cited from MathWorks*  

#### Anti-aliasing (Super-sampling AA vs. Random-sampling AA)
In Computer Graphics, aliasing refers to the visual artifacts that affect the quality of the images. In our particular case of ray-traced images, aliasing is presented as sharp step-like outlining of the 
rendered 3D objects. Aliasing occurs primarily because of insufficient sampling of the high-frequency details in the scene. The same concept lies behind Nyquist Theorem, a famous theorem in digital signal processing 
that states the relation between sample rate and signal bandwidth. The address aliasing, we have come up with 2 anti-aliasing (AA) algorithms. The first method is called super-sampling AA, in which each grid in front 
of the camera position is subdivided into many smaller grids, and a ray is fired through each subdivided grid. The final result of the original grid will be the average of the colors given by all the subdivided grids. 
This method is kind of a brute-force approach. It yields high-quality result but is more time-consuming. An alternative method we have come up with is the Random-sampling AA, in which we fire multiple rays through a 
grid such that their directions are distributed as a normal distribution. This method runs faster under equivalent hardware but produces more aliasing compared to super-sampling AA. We're currently still experimenting 
with these methods to try to find an equilibrium.  
![SSAA](readme_images/spheres_casting_soft_shadows_on_each_other_SSAA.png)  
*3D Scene Produced by SSAA*  
![Random AA](readme_images/spheres_casting_soft_shadows_on_each_other_RandomAA.png)  
*Same 3D Scene Produced by Random AA*  

#### Soft Shadow
