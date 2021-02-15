/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Adit Bhartia
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//Standard library stuff I included
#include <cmath>
#include <iostream>
#include <vector>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 1920
#define HEIGHT 1080
//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

struct Ray
{
	double Xo, Yo, Zo;
	double Xd, Yd, Zd; 
};

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

using namespace std;

double degreeToRadians(double degrees) {
	return (degrees / 180.0) * ((double) M_PI);
}

vector<double> pixelCenter(int i, int j) {
	double a = (double)WIDTH/HEIGHT;
	double tanTerm = tan(degreeToRadians(fov/2));
	double widthFactor = (2*a*tanTerm/WIDTH);
	double heightFactor = (2*tanTerm/HEIGHT); 
	double leftX = -1*a*tanTerm + widthFactor*(i);
	double rightX = leftX + widthFactor;
	double bottomY = -1*tanTerm + heightFactor*(j);
	double topY = bottomY + heightFactor;
	vector<double> center(3);
	center[0] = (leftX + rightX)/2;
	center[1] = (bottomY + topY)/2;
	center[2] = -1;
	return center;
}

double dotProduct(double* u, double* v) {
	return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

void crossProduct(double* u, double* v, double* result) {
	result[0] = u[1]*v[2] - u[2]*v[1];
	result[1] = u[2]*v[0] - u[0]*v[2];
	result[2] = u[0]*v[1] - u[1]*v[0];
}

double* normalize(double* input, double* result) {
        double scale = sqrt(pow(input[0], 2) + pow(input[1], 2) + pow(input[2], 2));
	result[0] = input[0]/scale;
	result[1] = input[1]/scale;
	result[2] = input[2]/scale;
	return result;
}

//specifically for vertices
double* vertexSubtract(Vertex left, Vertex right, double* result) {
	result[0] = left.position[0] - right.position[0];
	result[1] = left.position[1] - right.position[1];
	result[2] = left.position[2] - right.position[2];
	return result;
}

//specifically for vectors
vector<double> vectorSubtract(vector<double> left, vector<double> right) {
	vector<double> result(3);
	result[0] = left[0] - right[0];
	result[1] = left[1] - right[1];
	result[2] = left[2] - right[2];
	return result;
}

vector<double> evaluateRayAtPoint(Ray ray, double t) {
	vector<double> intersectionPoint{ray.Xo + ray.Xd*t, ray.Yo + ray.Yd*t, ray.Zo + ray.Zd*t};
	return intersectionPoint;
}

double xyProjectionArea(double Ax, double Ay, double Bx, double By, double Cx, double Cy) {
	return 0.5*((Bx-Ax)*(Cy-Ay) - (Cx-Ax)*(By-Ay));
}

double yzProjectionArea(double Ay, double Az, double By, double Bz, double Cy, double Cz) {
	return 0.5*((By-Ay)*(Cz-Az) - (Cy-Ay)*(Bz-Az));
}

double xzProjectionArea(double Ax, double Az, double Bx, double Bz, double Cx, double Cz) {
	return 0.5*((Bx-Ax)*(Cz-Az) - (Cx-Ax)*(Bz-Az));
}

void normalizeRayDirection(Ray& ray) {
        double scale = sqrt(pow(ray.Xd, 2) + pow(ray.Yd, 2) + pow(ray.Zd, 2));
	ray.Xd = ray.Xd/scale;
	ray.Yd = ray.Yd/scale;
	ray.Zd = ray.Zd/scale;
}

//using distance formula for seeing how close an intersection is
double rayToIntersectionDistance(Ray ray, vector<double> intersection) {
	int inside = pow(intersection[0] - ray.Xo,2) + pow(intersection[1] - ray.Yo,2) + pow(intersection[2] - ray.Zo,2);
	return sqrt(inside);
}

vector<double> calculateBarycentricCoords(Triangle triangle, vector<double> intersectionPoint) {
	Vertex a = triangle.v[0];
	Vertex b = triangle.v[1];
	Vertex c = triangle.v[2];
	double crossProductResult[3];
	double bMinusA[3];
	double cMinusA[3];
	crossProduct(vertexSubtract(b, a, bMinusA), vertexSubtract(c, a, cMinusA), crossProductResult);
	double normal[3];
	normalize(crossProductResult, normal);
	double xyPlane[] = {0, 0, 1};
	double xzPlane[] = {0, 1, 0};
	double yzPlane[] = {1, 0, 0};
	double denom = -1;
	double alpha = -1;
	double beta = -1;
	double gamma = -1;

	//find a plane that the normal is not perpendicular to and project to it to get area			
	if (dotProduct(normal, xzPlane) != 0) {  //projecting to xz plane
		denom = xzProjectionArea(a.position[0], a.position[2], b.position[0], b.position[2], c.position[0], c.position[2]);
		alpha = xzProjectionArea(intersectionPoint[0], intersectionPoint[2], b.position[0], b.position[2], c.position[0], c.position[2])/denom;
		beta = xzProjectionArea(a.position[0], a.position[2], intersectionPoint[0], intersectionPoint[2], c.position[0], c.position[2])/denom;
		gamma = xzProjectionArea(a.position[0], a.position[2], b.position[0], b.position[2], intersectionPoint[0], intersectionPoint[2])/denom;
	}

	else if (dotProduct(normal, xyPlane) != 0) { //projecting to xy plane
		denom = xyProjectionArea(a.position[0], a.position[1], b.position[0], b.position[1], c.position[0], c.position[1]);
		alpha = xyProjectionArea(intersectionPoint[0], intersectionPoint[1], b.position[0], b.position[1], c.position[0], c.position[1])/denom;
		beta = xyProjectionArea(a.position[0], a.position[1], intersectionPoint[0], intersectionPoint[1], c.position[0], c.position[1])/denom;
		gamma = xyProjectionArea(a.position[0], a.position[1], b.position[0], b.position[1], intersectionPoint[0], intersectionPoint[1])/denom;
	}

	else if (dotProduct(normal, yzPlane) != 0) { //projecting to yz plane
		denom = yzProjectionArea(a.position[1], a.position[2], b.position[1], b.position[2], c.position[1], c.position[2]);
		alpha = yzProjectionArea(intersectionPoint[1], intersectionPoint[2], b.position[1], b.position[2], c.position[1], c.position[2])/denom;
		beta = yzProjectionArea(a.position[1], a.position[2], intersectionPoint[1], intersectionPoint[2], c.position[1], c.position[2])/denom;
		gamma = yzProjectionArea(a.position[1], a.position[2], b.position[1], b.position[2], intersectionPoint[1], intersectionPoint[2])/denom;
	}

	vector<double> greeks{alpha, beta, gamma};
	return greeks;

}

//returns location of intersection
vector<double> CheckTriangleIntersection(Triangle &triangle, Ray &ray) {
	//Triangle vertices
	Vertex a = triangle.v[0];
	Vertex b = triangle.v[1];
	Vertex c = triangle.v[2];
	//calculate normal
	double crossProductResult[3];
	double bMinusA[3];
	double cMinusA[3];
	crossProduct(vertexSubtract(b, a, bMinusA), vertexSubtract(c, a, cMinusA), crossProductResult);
	double normal[3];
	normalize(crossProductResult, normal);
	//calculate t
	double negativeA[] = {-1*a.position[0], -1*a.position[1], -1*a.position[2]};
	double d = dotProduct(normal, negativeA);
	double rayOrigin[] = {ray.Xo, ray.Yo, ray.Zo};
	double rayDirection[] = {ray.Xd, ray.Yd, ray.Zd};
	double bottomTerm = dotProduct(normal, rayDirection);
	if (bottomTerm == 0) {
		return vector<double>();
	}
	double topTerm = -1.0*(dotProduct(normal, rayOrigin) + d);
	double t = topTerm/bottomTerm;
	if (t <= 0) return vector<double>();
	//Check if intersection point w/ plane inside triangle
	vector <double> intersectionPoint = evaluateRayAtPoint(ray, t);
	double xyPlane[] = {0, 0, 1};
	double xzPlane[] = {0, 1, 0};
	double yzPlane[] = {1, 0, 0};
	double denom, alpha, beta, gamma;

	//find a plane that normal is not perpendicular to, and project for calculating area
	if (dotProduct(normal, xzPlane) != 0) { 
		denom = xzProjectionArea(a.position[0], a.position[2], b.position[0], b.position[2], c.position[0], c.position[2]);
		alpha = xzProjectionArea(intersectionPoint[0], intersectionPoint[2], b.position[0], b.position[2], c.position[0], c.position[2])/denom;
		beta = xzProjectionArea(a.position[0], a.position[2], intersectionPoint[0], intersectionPoint[2], c.position[0], c.position[2])/denom;
		gamma = xzProjectionArea(a.position[0], a.position[2], b.position[0], b.position[2], intersectionPoint[0], intersectionPoint[2])/denom;
	}

	else if (dotProduct(normal, xyPlane) != 0) { 
		denom = xyProjectionArea(a.position[0], a.position[1], b.position[0], b.position[1], c.position[0], c.position[1]);
		alpha = xyProjectionArea(intersectionPoint[0], intersectionPoint[1], b.position[0], b.position[1], c.position[0], c.position[1])/denom;
		beta = xyProjectionArea(a.position[0], a.position[1], intersectionPoint[0], intersectionPoint[1], c.position[0], c.position[1])/denom;
		gamma = xyProjectionArea(a.position[0], a.position[1], b.position[0], b.position[1], intersectionPoint[0], intersectionPoint[1])/denom;
	}

	else if (dotProduct(normal, yzPlane) != 0) {
		denom = yzProjectionArea(a.position[1], a.position[2], b.position[1], b.position[2], c.position[1], c.position[2]);
		alpha = yzProjectionArea(intersectionPoint[1], intersectionPoint[2], b.position[1], b.position[2], c.position[1], c.position[2])/denom;
		beta = yzProjectionArea(a.position[1], a.position[2], intersectionPoint[1], intersectionPoint[2], c.position[1], c.position[2])/denom;
		gamma = yzProjectionArea(a.position[1], a.position[2], b.position[1], b.position[2], intersectionPoint[1], intersectionPoint[2])/denom;
	} 

	if (alpha < 0 || beta < 0 || gamma < 0) {
		return vector<double>();
	}
	if (alpha > 1 || beta > 1 || gamma > 1) {
		return vector<double>();
	}

	else return intersectionPoint;
}

vector<double> CheckSphereIntersection(Sphere sphere, Ray ray) {
	double Xo = ray.Xo;
	double Yo = ray.Yo;
	double Zo = ray.Zo;
	double Xd = ray.Xd;
	double Yd = ray.Yd;
	double Zd = ray.Zd;
	double Xc = sphere.position[0];
	double Yc = sphere.position[1];
	double Zc = sphere.position[2];
	double r = sphere.radius;
	
	double smallValue = 0.001;
	double a = 1.0;
	double b = 2*(Xd*(Xo-Xc) + Yd*(Yo-Yc) + Zd*(Zo-Zc));
	double c = pow(Xo-Xc, 2) + pow(Yo-Yc, 2) + pow(Zo-Zc, 2) - pow(r, 2);
	if ((pow(b, 2) -4*a*c) <= 0) {
		return vector<double>();
	}
	double t0 = (-1*b + sqrt(pow(b, 2) - 4*a*c))/2;
	double t1 = (-1*b - sqrt(pow(b, 2) - 4*a*c))/2;
	if (t0 > smallValue && t1 > smallValue) { //if both parameters are valid, return min
		return evaluateRayAtPoint(ray, min(t0, t1));
	}
	else if (t0 > smallValue) return evaluateRayAtPoint(ray, t0);
	else if (t1 > smallValue) return evaluateRayAtPoint(ray, t1);
	else return vector<double>();
}

//true if blocked by another object, false if ray to light from intersection is unblocked
bool shadowExists(Light light, vector<double> originalIntersection) {
	struct Ray shadowRay, n;
	shadowRay.Xo = originalIntersection[0];
	shadowRay.Yo = originalIntersection[1];
	shadowRay.Zo = originalIntersection[2];
	vector<double> lightPosition{light.position[0], light.position[1], light.position[2]};
	vector<double> rayDirection = vectorSubtract(lightPosition, originalIntersection);	
	shadowRay.Xd = rayDirection[0];
	shadowRay.Yd = rayDirection[1];
	shadowRay.Zd = rayDirection[2];
	normalizeRayDirection(shadowRay);
	//double distanceToLight = sqrt(pow(rayDirection[0], 2) + pow(rayDirection[1], 2) + pow(rayDirection[2], 2));
	double distanceToLight = pow(rayDirection[0], 2) + pow(rayDirection[1], 2) + pow(rayDirection[2], 2);
	double eps = 0.0001;
	//for all objects
	for (int i = 0; i<num_spheres; i++) {

		vector<double> intersectionPoint = CheckSphereIntersection(spheres[i], shadowRay);
		if (intersectionPoint.empty()) continue; //shadow ray does not hit this sphere
		double distanceToSphere = pow(intersectionPoint[0] - originalIntersection[0], 2) + pow(intersectionPoint[1] - originalIntersection[1], 2) + pow(intersectionPoint[2] - originalIntersection[2], 2);
		if (distanceToSphere > eps && distanceToSphere < distanceToLight) return true; //if sphere in the way of light, shadow exists
	}

	for (int i = 0; i<num_triangles; i++) {

		vector<double> intersectionPoint = CheckTriangleIntersection(triangles[i], shadowRay);
		if (intersectionPoint.empty()) continue; //shadow ray does not hit this triangle
		double distanceToTriangle = pow(intersectionPoint[0] - originalIntersection[0], 2) + pow(intersectionPoint[1] - originalIntersection[1], 2) + pow(intersectionPoint[2] - originalIntersection[2], 2);
		//potential further optimization: just store non-squre absolute version, since sqrt quite expensive to keep re-calculating
		if (distanceToTriangle > eps && distanceToTriangle < distanceToLight) return true; //if triangle in the way of light, shadow exists
	}
	return false;
}


//MODIFY THIS FUNCTION
void draw_scene()
{
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {
      vector<double> rayDirection = pixelCenter(x, y);
      struct Ray primaryRay;
      primaryRay.Xo = 0;
      primaryRay.Yo = 0;
      primaryRay.Zo = 0;
      primaryRay.Xd = rayDirection[0];
      primaryRay.Yd = rayDirection[1];
      primaryRay.Zd = rayDirection[2];
      normalizeRayDirection(primaryRay);

      double color[] = {ambient_light[0], ambient_light[1], ambient_light[2]};
      //variables to keep track of which is the closest object so we can apply shadow ray + phong for only that closest intersecting object
      double closestDistance = INT_MAX; 
      Sphere closestSphere;
      Triangle closestTriangle;
      bool closestIsSphere;
      vector<double> closestIntersection;
      //if not set to true, means that pixel's primary ray did not hit any objects and is the background
      bool objectFound = false;
	for (int i = 0; i<num_spheres; i++) {
		vector<double> SphereIntersection = CheckSphereIntersection(spheres[i], primaryRay);
		if (SphereIntersection.empty()) continue;
		else {
			double distanceToSphere = rayToIntersectionDistance(primaryRay, SphereIntersection);
			if (distanceToSphere < closestDistance) { //update vars if found a closer object
				closestDistance = distanceToSphere;
				closestIntersection = SphereIntersection;
				closestSphere = spheres[i];
				closestIsSphere = true;		
				objectFound = true;
			}
		}
	}

	for (int i = 0; i<num_triangles; i++) {
		vector<double> TriangleIntersection = CheckTriangleIntersection(triangles[i], primaryRay);
		if (TriangleIntersection.empty()) continue;
		else {
			double distanceToTriangle = rayToIntersectionDistance(primaryRay, TriangleIntersection);
			if (distanceToTriangle < closestDistance) { //update vars if found a closer object
				closestDistance = distanceToTriangle;
				closestIntersection = TriangleIntersection;
				closestTriangle = triangles[i];
				closestIsSphere = false;		
				objectFound = true;
			}
		}
	}

	if (objectFound) { //if triangle/sphere found with primary ray
		for (int i = 0; i<num_lights; i++) {
			if (!shadowExists(lights[i], closestIntersection)) {
				//Common calculations for both sphere + triangle phong shading
				double vecToLight[3];
				vecToLight[0] = lights[i].position[0] - closestIntersection[0];
				vecToLight[1] = lights[i].position[1] - closestIntersection[1];
				vecToLight[2] = lights[i].position[2] - closestIntersection[2];
				double l[3];
				normalize(vecToLight, l);
				double v[3];//vector to camera (primary ray origin) 
				v[0] = primaryRay.Xo - closestIntersection[0];
				v[1] = primaryRay.Yo - closestIntersection[1];
				v[2] = primaryRay.Zo - closestIntersection[2];
				normalize(v, v);	
				double localColor[3];
				//if the closest intersection is a sphere
				if (closestIsSphere) {
					//Phong for sphere
					double sphereNormal[3];
					sphereNormal[0] = closestIntersection[0] - closestSphere.position[0];
					sphereNormal[1] = closestIntersection[1] - closestSphere.position[1];
					sphereNormal[2] = closestIntersection[2] - closestSphere.position[2];
					double n[3];
					normalize(sphereNormal, n);
					double lDotN = dotProduct(l, n);
					if (lDotN < 0) lDotN = 0;
					//reflected ray
					double r[3];
					r[0] = 2*lDotN*n[0] - l[0];
					r[1] = 2*lDotN*n[1] - l[1];
					r[2] = 2*lDotN*n[2] - l[2];
					normalize(r, r);
					double rDotV = dotProduct(r, v);
					if (rDotV < 0) rDotV = 0;
					//color for each channel
					localColor[0] = lights[i].color[0]*(closestSphere.color_diffuse[0]*(lDotN) + closestSphere.color_specular[0]*(pow(rDotV, closestSphere.shininess)));
					localColor[1] = lights[i].color[1]*(closestSphere.color_diffuse[1]*(lDotN) + closestSphere.color_specular[1]*(pow(rDotV, closestSphere.shininess)));
					localColor[2] = lights[i].color[2]*(closestSphere.color_diffuse[2]*(lDotN) + closestSphere.color_specular[2]*(pow(rDotV, closestSphere.shininess)));

				}

				else { //closest is triangle
					//phong for triangle
					vector<double> greeks = calculateBarycentricCoords(closestTriangle, closestIntersection);
					double alpha = greeks[0];
					double beta = greeks[1];
					double gamma = greeks[2];
					Vertex a = closestTriangle.v[0];
					Vertex b = closestTriangle.v[1];
					Vertex c = closestTriangle.v[2];
					//interpolate coordinates of normal, and then normalize
					//interpolate diffuse, specular, and shininess
					double intersectionNormal[3];
					double intersectionDiffuse[3];
					double intersectionSpecular[3];
					//Interpolated values for given intersection
					for (int colorIndex = 0; colorIndex<3; colorIndex++) {
						intersectionNormal[colorIndex] = alpha*a.normal[colorIndex] + beta*b.normal[colorIndex] + gamma*c.normal[colorIndex];
						intersectionDiffuse[colorIndex] = alpha*a.color_diffuse[colorIndex] + beta*b.color_diffuse[colorIndex] + gamma*c.color_diffuse[colorIndex];
						intersectionSpecular[colorIndex] = alpha*a.color_specular[colorIndex] + beta*b.color_specular[colorIndex] + gamma*c.color_specular[colorIndex];
					}

					double shininess = alpha*a.shininess + beta*b.shininess + gamma*c.shininess;
					double n[3];
					normalize(intersectionNormal, n);
					double lDotN = dotProduct(l, n);
					if (lDotN < 0) lDotN = 0;
					//reflected ray
					double r[3];
					r[0] = 2*lDotN*n[0] - l[0];
					r[1] = 2*lDotN*n[1] - l[1];
					r[2] = 2*lDotN*n[2] - l[2];
					normalize(r, r);
					double rDotV = dotProduct(r, v);
					if (rDotV < 0) rDotV = 0;
					//color for each channel
					localColor[0] = lights[i].color[0]*(intersectionDiffuse[0]*(lDotN) + intersectionSpecular[0]*(pow(rDotV, shininess)));
					localColor[1] = lights[i].color[1]*(intersectionDiffuse[1]*(lDotN) + intersectionSpecular[1]*(pow(rDotV, shininess)));
					localColor[2] = lights[i].color[2]*(intersectionDiffuse[2]*(lDotN) + intersectionSpecular[2]*(pow(rDotV, shininess)));	
				}
				//update color of pixel if not in shadow w/ respect to that light
				color[0] += localColor[0];
				color[1] += localColor[1];
				color[2] += localColor[2];
			} 
		}
	}
	//clamping
	if (color[0] > 1) color[0] = 1;
	if (color[1] > 1) color[1] = 1;
	if (color[2] > 1) color[2] = 1;
	//background color is white
      if (!objectFound) plot_pixel(x, y, 255, 255, 255);
      else plot_pixel(x, y, round(color[0]*255.0), round(color[1]*255.0), round(color[2]*255.0));
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);
    
      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

