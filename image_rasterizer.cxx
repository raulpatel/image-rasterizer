// Access my movie at: https://youtu.be/aKMaKMug-ZU
#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <cmath>

using std::cerr;
using std::endl;

/********************
 * Helper functions *
 ********************/
double ceiling_function(double f) {
    return ceil(f - 0.00001);
}

double floor_function(double f) {
    return floor(f + 0.00001);
}

void CrossProduct(double * A, double * B, double( & prod)[3]) {
    // helper function for cross product
    prod[0] = (A[1] * B[2]) - (A[2] * B[1]);
    prod[1] = (A[2] * B[0]) - (A[0] * B[2]);
    prod[2] = (A[0] * B[1]) - (A[1] * B[0]);
}

double DotProduct(double * A, double * B) {
    // helper function for dot product
    return (A[0] * B[0]) + (A[1] * B[1]) + (A[2] * B[2]);
}

double FindNorm(double * A) {
    // helper function for finding the norm of a vector
    return sqrt((A[0] * A[0]) + (A[1] * A[1]) + (A[2] * A[2]));
}

void NormalizeVect(double( & A)[3]) {
    // helper function for normalizing vectors
    double norm = FindNorm(A);
    if (norm > 0) {
        for (int i = 0; i < 3; i++) {
            A[i] /= norm;
        }
    }
}

/*****************
 * VTK functions *
 *****************/
vtkImageData *
NewImage(int height, int width) {
    vtkImageData * img = vtkImageData::New();
    img -> SetDimensions(width, height, 1);
    img -> AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData * img,
    const char * filename) {
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter * writer = vtkPNGWriter::New();
    writer -> SetInputData(img);
    writer -> SetFileName(full_filename.c_str());
    writer -> Write();
    writer -> Delete();
}

/*****************************
 * Class/struct declarations *
 *****************************/
class Screen {
    public:
        unsigned char * buffer;
        double * zbuffer;
        int width, height;
        void colorPixel(int idx, double r, double g,
                        double b, double s, double z);
};

struct LightingParameters {
    LightingParameters(void) {
        lightDir[0] = -0.6;
        lightDir[1] = 0;
        lightDir[2] = -0.8;
        Ka = 0.3;
        Kd = 0.7;
        Ks = 2.8;
        alpha = 50.5;
    };

    double lightDir[3]; // The direction of the light source
    double Ka; // The coefficient for ambient lighting
    double Kd; // The coefficient for diffuse lighting
    double Ks; // The coefficient for specular lighting
    double alpha; // The exponent term for specular lighting
};

class Camera;

class Triangle {
    public:
        double X[3];
    double Y[3];
    double Z[3];
    double colors[3][3];
    double normals[3][3];
    double shaders[3];
    // markers for each vertex
    int left, mid, right;
    // helper function to split an arb tri into two, pass left and right tris by reference
    void splitTriangle(Triangle & lt, Triangle & rt);
    // function to rasterize an arb tri by splitting and rasterizing the sub tris
    void rasterizeArbitraryTriangle(Screen screen);
    // helper functions for each type of triangle to rasterize, called by rastArb
    void rasterizeLeftTriangle(Screen screen);
    void rasterizeRightTriangle(Screen screen);
    // helper function to calculate the shading for each point in the triangle
    void calculateShading(LightingParameters & lp, Camera & c);
};

class Matrix {
    public:
        double A[4][4]; // A[i][j] means row i, column j

    void TransformPoint(const double * ptIn, double * ptOut);
    static Matrix ComposeMatrices(const Matrix & ,
        const Matrix & );
    void Print(ostream & o);
};

typedef struct CameraFrame {
    double u[3];
    double v[3];
    double w[3];
    double O[3];
}
CameraFrame;

class Camera {
    public:
        double near, far;
    double angle;
    double position[3];
    double focus[3];
    double up[3];
    CameraFrame frame;

    Matrix ViewTransform(void);
    Matrix CameraTransform(void);
    Matrix DeviceTransform(void);
// matrix multiplications to transform points from world space to device space
    std::vector < Triangle > TransformTrisToDevSpace(
        std::vector < Triangle > tris);
};

/************************
 * Matrix Class methods *
 ************************/

void
Matrix::Print(ostream & o) {
    for (int i = 0; i < 4; i++) {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix
Matrix::ComposeMatrices(const Matrix & M1,
    const Matrix & M2) {
    Matrix rv;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            rv.A[i][j] = 0;
            for (int k = 0; k < 4; k++)
                rv.A[i][j] += M1.A[i][k] * M2.A[k][j];
        }

    return rv;
}

void
Matrix::TransformPoint(const double * ptIn, double * ptOut) {
    ptOut[0] = ptIn[0] * A[0][0] +
        ptIn[1] * A[1][0] +
        ptIn[2] * A[2][0] +
        ptIn[3] * A[3][0];
    ptOut[1] = ptIn[0] * A[0][1] +
        ptIn[1] * A[1][1] +
        ptIn[2] * A[2][1] +
        ptIn[3] * A[3][1];
    ptOut[2] = ptIn[0] * A[0][2] +
        ptIn[1] * A[1][2] +
        ptIn[2] * A[2][2] +
        ptIn[3] * A[3][2];
    ptOut[3] = ptIn[0] * A[0][3] +
        ptIn[1] * A[1][3] +
        ptIn[2] * A[2][3] +
        ptIn[3] * A[3][3];
}

/************************
 * Camera Class methods *
 ************************/

Matrix Camera::CameraTransform() {
    // set Camera Frame
    for (int i = 0; i < 3; i++) {
        frame.O[i] = position[i];
        frame.w[i] = position[i] - focus[i];
    }
    CrossProduct(up, frame.w, frame.u);
    // we have to normalize u, v, and w before use
    NormalizeVect(frame.u);
    CrossProduct(frame.w, frame.u, frame.v);
    NormalizeVect(frame.v);
    NormalizeVect(frame.w);

    Matrix ret;
    // initialize matrix to 0s
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ret.A[i][j] = 0;
        }
    }
    // set matrix values
    double t[3];
    for (int i = 0; i < 3; i++) {
        ret.A[i][0] = frame.u[i];
        ret.A[i][1] = frame.v[i];
        ret.A[i][2] = frame.w[i];
        ret.A[i][3] = 0;
        t[i] = 0 - frame.O[i];
    }
    ret.A[3][0] = DotProduct(frame.u, t);
    ret.A[3][1] = DotProduct(frame.v, t);
    ret.A[3][2] = DotProduct(frame.w, t);
    ret.A[3][3] = 1;
    return ret;
}

Matrix Camera::ViewTransform() {
    Matrix ret;
    // initialize matrix to 0s
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ret.A[i][j] = 0;
        }
    }
    // set matrix
    ret.A[0][0] = 1 / tan(angle / 2); // (1 / tan) = cot
    ret.A[1][1] = 1 / tan(angle / 2);
    ret.A[2][2] = (far + near) / (far - near);
    ret.A[2][3] = -1;
    ret.A[3][2] = 2 * far * near / (far - near);

    return ret;
}

Matrix Camera::DeviceTransform() {
    Matrix ret;
    // initialize matrix to 0s
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ret.A[i][j] = 0;
        }
    }
    // set matrix
    ret.A[0][0] = 500;
    ret.A[1][1] = 500;
    ret.A[2][2] = 1;
    ret.A[3][0] = 500;
    ret.A[3][1] = 500;
    ret.A[3][3] = 1;

    return ret;
}

double SineParameterize(int curFrame, int nFrames, int ramp) {
    int nNonRamp = nFrames - 2 * ramp;
    double height = 1. / (nNonRamp + 4 * ramp / M_PI);
    if (curFrame < ramp) {
        double factor = 2 * height * ramp / M_PI;
        double eval = cos(M_PI / 2 * ((double) curFrame) / ramp);
        return (1. - eval) * factor;
    } else if (curFrame > nFrames - ramp) {
        int amount_left = nFrames - curFrame;
        double factor = 2 * height * ramp / M_PI;
        double eval = cos(M_PI / 2 * ((double) amount_left / ramp));
        return 1. - (1 - eval) * factor;
    }
    double amount_in_quad = ((double) curFrame - ramp);
    double quad_part = amount_in_quad * height;
    double curve_part = height * (2 * ramp) / M_PI;
    return quad_part + curve_part;
}

Camera
GetCamera(int frame, int nframes) {
    double t = SineParameterize(frame, nframes, nframes / 10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI / 6;
    c.position[0] = 40 * sin(2 * M_PI * t);
    c.position[1] = 40 * cos(2 * M_PI * t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}

std::vector < Triangle >
    Camera::TransformTrisToDevSpace(std::vector < Triangle > tris) {
        std::vector < Triangle > newtris(tris.size());
        Matrix cam, view, dev, M;
        // camera transform
        cam = CameraTransform();
        // view transform
        view = ViewTransform();
        // device transform
        dev = DeviceTransform();
        // compose matrices so that we get one matrix that correctly sets points
        M = Matrix::ComposeMatrices(cam, view);
        M = Matrix::ComposeMatrices(M, dev);

        for (int i = 0; i < tris.size(); i++) {

            int left, mid, right;
            double ileft[4], imid[4], iright[4], oleft[4], omid[4], oright[4];
            left = tris[i].left;
            mid = tris[i].mid;
            right = tris[i].right;
            // set up arrays for use in TransformPoints
            ileft[0] = tris[i].X[left];
            ileft[1] = tris[i].Y[left];
            ileft[2] = tris[i].Z[left];
            ileft[3] = 1;
            imid[0] = tris[i].X[mid];
            imid[1] = tris[i].Y[mid];
            imid[2] = tris[i].Z[mid];
            imid[3] = 1;
            iright[0] = tris[i].X[right];
            iright[1] = tris[i].Y[right];
            iright[2] = tris[i].Z[right];
            iright[3] = 1;
            // transform points using the composed matrix
            M.TransformPoint(ileft, oleft);
            M.TransformPoint(imid, omid);
            M.TransformPoint(iright, oright);
            // divide all values by w to get the correct vals
            for (int j = 0; j < 4; j++) {
                oleft[j] /= oleft[3];
                omid[j] /= omid[3];
                oright[j] /= oright[3];
            }
            // set the new points
            newtris[i].X[left] = oleft[0];
            newtris[i].X[mid] = omid[0];
            newtris[i].X[right] = oright[0];
            newtris[i].Y[left] = oleft[1];
            newtris[i].Y[mid] = omid[1];
            newtris[i].Y[right] = oright[1];
            newtris[i].Z[left] = oleft[2];
            newtris[i].Z[mid] = omid[2];
            newtris[i].Z[right] = oright[2];
            // set shader values before converting vertices
            newtris[i].shaders[0] = tris[i].shaders[0];
            newtris[i].shaders[1] = tris[i].shaders[1];
            newtris[i].shaders[2] = tris[i].shaders[2];

            // adjust left, mid, and right markers
            newtris[i].left = std::min_element(newtris[i].X, newtris[i].X + 3) -
                newtris[i].X;
            newtris[i].right = std::max_element(newtris[i].X, newtris[i].X + 3) -
                newtris[i].X;
            newtris[i].mid = 3 - newtris[i].left - newtris[i].right;
            // set the colors
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    newtris[i].colors[j][k] = tris[i].colors[j][k];
                    newtris[i].normals[j][k] = tris[i].normals[j][k];
                }
            }
        }
        return newtris;
    }

/****************************
 * Get LightingParameters  *
 ****************************/

LightingParameters
GetLighting(Camera c) {
    LightingParameters lp;
    lp.lightDir[0] = c.position[0] - c.focus[0];
    lp.lightDir[1] = c.position[1] - c.focus[1];
    lp.lightDir[2] = c.position[2] - c.focus[2];
    double mag = sqrt(lp.lightDir[0] * lp.lightDir[0] +
        lp.lightDir[1] * lp.lightDir[1] +
        lp.lightDir[2] * lp.lightDir[2]);
    if (mag > 0) {
        lp.lightDir[0] /= mag;
        lp.lightDir[1] /= mag;
        lp.lightDir[2] /= mag;
    }

    return lp;
}

/***********************
 * Screen Class method *
 ***********************/

void Screen::colorPixel(int idx, double r, double g, double b, double s, double z) {
    // make sure that the pixel is in frame and in front
    if (idx >= 0 && idx < 1000000 && z > zbuffer[idx]) {
        int curr = idx * 3;
        buffer[curr] = ceiling_function(255 * std::min(1.0, (r * s)));
        buffer[curr + 1] = ceiling_function(255 * std::min(1.0, (g * s)));
        buffer[curr + 2] = ceiling_function(255 * std::min(1.0, (b * s)));
        zbuffer[idx] = z;
    }
}

/**************************
 * Triangle Class methods *
 **************************/

void Triangle::calculateShading(LightingParameters & lp, Camera & c) {
    // shading amount = ka + kd * diffuse + ks * specular
    // = (ka + kd * (L . N) + ks * (V . (2 * (L . N) * N - L))^alpha)
    // (L . N) dot product of light direction and normal of the vertex
    double ldot0 = DotProduct(lp.lightDir, normals[0]);
    double ldot1 = DotProduct(lp.lightDir, normals[1]);
    double ldot2 = DotProduct(lp.lightDir, normals[2]);
    // reflection value AKA 2 * (L . N) * N - L
    double r0[3] = {
        (2 * ldot0 * normals[0][0] - lp.lightDir[0]),
        (2 * ldot0 * normals[0][1] - lp.lightDir[1]),
        (2 * ldot0 * normals[0][2] - lp.lightDir[2])
    };
    double r1[3] = {
        (2 * ldot1 * normals[1][0] - lp.lightDir[0]),
        (2 * ldot1 * normals[1][1] - lp.lightDir[1]),
        (2 * ldot1 * normals[1][2] - lp.lightDir[2])
    };
    double r2[3] = {
        (2 * ldot2 * normals[2][0] - lp.lightDir[0]),
        (2 * ldot2 * normals[2][1] - lp.lightDir[1]),
        (2 * ldot2 * normals[2][2] - lp.lightDir[2])
    };
    // now check that (L . N) >= 0
    ldot0 = std::max(0.0, ldot0);
    ldot1 = std::max(0.0, ldot1);
    ldot2 = std::max(0.0, ldot2);
    // find view direction
    double v0[3] = {
        (c.position[0] - X[0]),
        (c.position[1] - Y[0]),
        (c.position[2] - Z[0])
    };
    double v1[3] = {
        (c.position[0] - X[1]),
        (c.position[1] - Y[1]),
        (c.position[2] - Z[1])
    };
    double v2[3] = {
        (c.position[0] - X[2]),
        (c.position[1] - Y[2]),
        (c.position[2] - Z[2])
    };
    // Normalize vectors
    NormalizeVect(v0);
    NormalizeVect(v1);
    NormalizeVect(v2);
    NormalizeVect(r0);
    NormalizeVect(r1);
    NormalizeVect(r2);
    // find specular lighting
    double vdot0 = pow(DotProduct(r0, v0), lp.alpha);
    double vdot1 = pow(DotProduct(r1, v1), lp.alpha);
    double vdot2 = pow(DotProduct(r2, v2), lp.alpha);
    // if pow function returns NaN change to 0
    if (std::isnan(vdot0)) vdot0 = 0;
    if (std::isnan(vdot1)) vdot1 = 0;
    if (std::isnan(vdot2)) vdot2 = 0;
    // set shaders
    shaders[0] = lp.Ka + lp.Kd * ldot0 + lp.Ks * vdot0;
    shaders[1] = lp.Ka + lp.Kd * ldot1 + lp.Ks * vdot1;
    shaders[2] = lp.Ka + lp.Kd * ldot2 + lp.Ks * vdot2;

}

void Triangle::splitTriangle(Triangle & lt, Triangle & rt) {
    // find slope of edge opposite middle vertex 
    double bot_slope = (Y[right] - Y[left]) / (X[right] - X[left]);
    double t = (X[mid] - X[left]) / (X[right] - X[left]);
    // use the slope to find the y value of our new third vertex
    double center_y = Y[left] + (bot_slope * (X[mid] - X[left]));
    // LERP (linear interpolation to estimate z value and color)
    double center_z = Z[left] + (t * (Z[right] - Z[left]));
    double center_r = colors[left][0] + (t * (colors[right][0] - colors[left][0]));
    double center_g = colors[left][1] + (t * (colors[right][1] - colors[left][1]));
    double center_b = colors[left][2] + (t * (colors[right][2] - colors[left][2]));
    double center_shader = shaders[left] + (t * (shaders[right] - shaders[left]));
    // initialize coordinates and color for left triangle
    lt.X[left] = X[left];
    lt.X[mid] = X[mid];
    lt.X[right] = X[mid];

    lt.Y[left] = Y[left];
    lt.Y[mid] = Y[mid];
    lt.Y[right] = center_y;

    lt.Z[left] = Z[left];
    lt.Z[mid] = Z[mid];
    lt.Z[right] = center_z;

    lt.shaders[left] = shaders[left];
    lt.shaders[mid] = shaders[mid];
    lt.shaders[right] = center_shader;

    lt.left = left;
    lt.mid = mid;
    lt.right = right;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            lt.colors[i][j] = colors[i][j];
        }
        if (i == right) {
            lt.colors[i][0] = center_r;
            lt.colors[i][1] = center_g;
            lt.colors[i][2] = center_b;
        }
    }

    // same but right triangle
    rt.X[left] = X[mid];
    rt.X[mid] = X[mid];
    rt.X[right] = X[right];

    rt.Y[left] = center_y;
    rt.Y[mid] = Y[mid];
    rt.Y[right] = Y[right];

    rt.Z[left] = center_z;
    rt.Z[mid] = Z[mid];
    rt.Z[right] = Z[right];

    rt.shaders[left] = center_shader;
    rt.shaders[mid] = shaders[mid];
    rt.shaders[right] = shaders[right];

    rt.left = left;
    rt.mid = mid;
    rt.right = right;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            rt.colors[i][j] = colors[i][j];
        }
        if (i == left) {
            rt.colors[i][0] = center_r;
            rt.colors[i][1] = center_g;
            rt.colors[i][2] = center_b;
        }
    }
}

void Triangle::rasterizeLeftTriangle(Screen screen) {

    double colLeft, colRight, top, bottom, ztop, zbot, rtop, rbot, gtop, gbot, btop, bbot, stop, sbot;

    colRight = X[right];
    int columnMax = floor_function(colRight);

    colLeft = X[left];
    int columnMin = ceiling_function(colLeft);
    // find the slope for the top and bottom for boundaries of inner loop
    double top_slope = (Y[mid] - Y[left]) / (X[mid] - X[left]);
    double bot_slope = (Y[right] - Y[left]) / (X[right] - X[left]);

    for (int i = columnMin; i <= columnMax; i++) {
        // check boundaries
        if (i < 0) continue;
        if (i > 999) break;
        double t = (i - X[left]) / (X[right] - X[left]); // ratio of field
        double inter1 = Y[left] + ((i - colLeft) * bot_slope); // intersection points
        double inter2 = Y[left] + ((i - colLeft) * top_slope);
        // LERP values for z,r,g,b at both endpoints
        double z1 = Z[left] + (t * (Z[right] - Z[left]));
        double z2 = Z[left] + (t * (Z[mid] - Z[left]));
        double s1 = shaders[left] + (t * (shaders[right] - shaders[left]));
        double s2 = shaders[left] + (t * (shaders[mid] - shaders[left]));
        double r1 = colors[left][0] + (t * (colors[right][0] - colors[left][0]));
        double r2 = colors[left][0] + (t * (colors[mid][0] - colors[left][0]));
        double g1 = colors[left][1] + (t * (colors[right][1] - colors[left][1]));
        double g2 = colors[left][1] + (t * (colors[mid][1] - colors[left][1]));
        double b1 = colors[left][2] + (t * (colors[right][2] - colors[left][2]));
        double b2 = colors[left][2] + (t * (colors[mid][2] - colors[left][2]));
        // check orientation before using
        if (inter1 > inter2) {
            top = inter1;
            ztop = z1;
            stop = s1;
            rtop = r1;
            gtop = g1;
            btop = b1;
            bottom = inter2;
            zbot = z2;
            sbot = s2;
            rbot = r2;
            gbot = g2;
            bbot = b2;
        } else {
            top = inter2;
            ztop = z2;
            stop = s2;
            rtop = r2;
            gtop = g2;
            btop = b2;
            bottom = inter1;
            zbot = z1;
            sbot = s1;
            rbot = r1;
            gbot = g1;
            bbot = b1;
        }
        int bottomEnd = ceiling_function(bottom);
        int topEnd = floor_function(top);

        for (int j = bottomEnd; j <= topEnd; j++) {
            double tvert;
            if (bottom == top) tvert = 0; // check for divide by 0
            else tvert = (j - bottom) / (top - bottom);
            // LERP again
            double zpix = zbot + (tvert * (ztop - zbot));
            double spix = sbot + (tvert * (stop - sbot));
            double red = rbot + (tvert * (rtop - rbot));
            double green = gbot + (tvert * (gtop - gbot));
            double blue = bbot + (tvert * (btop - bbot));
            int curr = (j * screen.width) + i;
            screen.colorPixel(curr, red, green, blue, spix, zpix);
        }
    }
}

void Triangle::rasterizeRightTriangle(Screen screen) {

    double colLeft, colRight, top, bottom, ztop, zbot, rtop, rbot, gtop, gbot, btop, bbot, stop, sbot;

    colRight = X[right];
    int columnMax = floor_function(colRight);

    colLeft = X[left];
    int columnMin = ceiling_function(colLeft);
    // find the slope for the top and bottom for boundaries of inner loop
    double top_slope = (Y[right] - Y[mid]) / (X[right] - X[mid]);
    double bot_slope = (Y[right] - Y[left]) / (X[right] - X[left]);

    for (int i = columnMin; i <= columnMax; i++) {
        // set boundaries
        if (i < 0) continue;
        if (i > 999) break;
        double t = (i - X[left]) / (X[right] - X[left]);
        double inter1 = Y[left] + ((i - colLeft) * bot_slope);
        double inter2 = Y[mid] + ((i - colLeft) * top_slope);
        // LERP for depth and color at both ends
        double z1 = Z[left] + (t * (Z[right] - Z[left]));
        double z2 = Z[mid] + (t * (Z[right] - Z[mid]));
        double s1 = shaders[left] + (t * (shaders[right] - shaders[left]));
        double s2 = shaders[mid] + (t * (shaders[right] - shaders[mid]));
        double r1 = colors[left][0] + (t * (colors[right][0] - colors[left][0]));
        double r2 = colors[mid][0] + (t * (colors[right][0] - colors[mid][0]));
        double g1 = colors[left][1] + (t * (colors[right][1] - colors[left][1]));
        double g2 = colors[mid][1] + (t * (colors[right][1] - colors[mid][1]));
        double b1 = colors[left][2] + (t * (colors[right][2] - colors[left][2]));
        double b2 = colors[mid][2] + (t * (colors[right][2] - colors[mid][2]));
        // check orientation before using
        if (inter1 > inter2) {
            top = inter1;
            ztop = z1;
            stop = s1;
            rtop = r1;
            gtop = g1;
            btop = b1;
            bottom = inter2;
            zbot = z2;
            sbot = s2;
            rbot = r2;
            gbot = g2;
            bbot = b2;
        } else {
            top = inter2;
            ztop = z2;
            stop = s2;
            rtop = r2;
            gtop = g2;
            btop = b2;
            bottom = inter1;
            zbot = z1;
            sbot = s1;
            rbot = r1;
            gbot = g1;
            bbot = b1;
        }
        int bottomEnd = ceiling_function(bottom);
        int topEnd = floor_function(top);

        for (int j = bottomEnd; j <= topEnd; j++) {
            double tvert;
            if (bottom == top) tvert = 0; // avoid divide by 0
            else tvert = (j - bottom) / (top - bottom);
            double zpix = zbot + (tvert * (ztop - zbot));
            double spix = sbot + (tvert * (stop - sbot));
            double red = rbot + (tvert * (rtop - rbot));
            double green = gbot + (tvert * (gtop - gbot));
            double blue = bbot + (tvert * (btop - bbot));
            int curr = (j * screen.width) + i;
            screen.colorPixel(curr, red, green, blue, spix, zpix);
        }

    }
}

void Triangle::rasterizeArbitraryTriangle(Screen screen) {
    Triangle lt, rt;
    splitTriangle(lt, rt);
    lt.rasterizeLeftTriangle(screen);
    rt.rasterizeRightTriangle(screen);
}

/***********************
 * Parse Triangle Data *
 ***********************/

std::vector < Triangle >
GetTriangleData(void) {
    vtkPolyDataReader * rdr = vtkPolyDataReader::New();
    rdr -> SetFileName("aneurysm_model.vtk");
    cerr << "Reading" << endl;
    rdr -> Update();
    cerr << "Done reading" << endl;
    if (rdr -> GetOutput() -> GetNumberOfCells() == 0) {
       cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData * pd = rdr -> GetOutput();

    int numTris = pd -> GetNumberOfCells();
    vtkPoints * pts = pd -> GetPoints();
    vtkCellArray * cells = pd -> GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray * ) pd -> GetPointData() -> GetArray("hardyglobal");
    double * color_ptr = var -> GetPointer(0);
        
    vtkFloatArray * n = (vtkFloatArray * ) pd -> GetPointData() -> GetNormals();
    float * normals = n -> GetPointer(0);
    std::vector < Triangle > tris(numTris);
    vtkIdType npts;
    vtkIdType * ptIds;
    int idx;
    for (idx = 0, cells -> InitTraversal(); cells -> GetNextCell(npts, ptIds); idx++) {
        if (npts != 3) {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double * pt = NULL;
        pt = pts -> GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
        
        tris[idx].normals[0][0] = normals[3 * ptIds[0] + 0];
        tris[idx].normals[0][1] = normals[3 * ptIds[0] + 1];
        tris[idx].normals[0][2] = normals[3 * ptIds[0] + 2];
        
        pt = pts -> GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
        
        tris[idx].normals[1][0] = normals[3 * ptIds[1] + 0];
        tris[idx].normals[1][1] = normals[3 * ptIds[1] + 1];
        tris[idx].normals[1][2] = normals[3 * ptIds[1] + 2];
        
        pt = pts -> GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
        
        tris[idx].normals[2][0] = normals[3 * ptIds[2] + 0];
        tris[idx].normals[2][1] = normals[3 * ptIds[2] + 1];
        tris[idx].normals[2][2] = normals[3 * ptIds[2] + 2];
        
        // analyze vertices
        tris[idx].left = std::min_element(tris[idx].X, tris[idx].X + 3) -
            tris[idx].X;
        tris[idx].right = std::max_element(tris[idx].X, tris[idx].X + 3) -
            tris[idx].X;
        tris[idx].mid = 3 - tris[idx].left - tris[idx].right;

        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 },
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 } 
                                  };
        for (int j = 0; j < 3; j++) {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0; r < 7; r++) {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7) {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val - mins[r]) / (maxs[r] - mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0] + proportion * (RGB[r + 1][0] - RGB[r][0])) / 255.0;
            tris[idx].colors[j][1] = (RGB[r][1] + proportion * (RGB[r + 1][1] - RGB[r][1])) / 255.0;
            tris[idx].colors[j][2] = (RGB[r][2] + proportion * (RGB[r + 1][2] - RGB[r][2])) / 255.0;
        }
    }

    return tris;
}

/********
 * Main *
 ********/

int main() {
    vtkImageData * image = NewImage(1000, 1000);
    unsigned char * buffer =
    (unsigned char * ) image -> GetScalarPointer(0, 0, 0);
    int npixels = 1000 * 1000;
    double zbuffer[npixels];
    std::vector < Triangle > triangles = GetTriangleData();

    for (int i = 0; i < 1 /*1000*/ ; i++) { // set to 1 for submission
        int curr_frame = i;

        for (int i = 0; i < npixels; i++)
            zbuffer[i] = -1;
        for (int i = 0; i < npixels * 3; i++)
            buffer[i] = 0;

        Screen screen;
        screen.buffer = buffer;
        screen.zbuffer = zbuffer;
        screen.width = 1000;
        screen.height = 1000;

        // initialize camera
        Camera c = GetCamera(curr_frame, 1000);
        LightingParameters lp = GetLighting(c);
        // find shading before altering vertices
        for (int j = 0; j < triangles.size(); j++) {
            triangles[j].calculateShading(lp, c);
        }
        // adjust the points relative to the angle
        std::vector < Triangle > tris = c.TransformTrisToDevSpace(triangles);
        //NormalizeVect(c.position);

        // YOUR CODE GOES HERE TO DEPOSIT TRIANGLES INTO PIXELS USING THE SCANLINE ALGORITHM
        for (int j = 0; j < tris.size(); j++) {
            tris[j].rasterizeArbitraryTriangle(screen);
        }
        char name[16];
        sprintf(name, "frames/frame%03d", i);
        WriteImage(image, name);
    }
}
