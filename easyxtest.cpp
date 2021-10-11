// easyxtest.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#define _USE_MATH_DEFINES
#define _AFXDLL
#define  _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <graphics.h>
#include <time.h>
#include<math.h>
#include<iostream>
#include<cstring>
#include<vector>
#include<queue>
#include<string>
#include<algorithm>
//#include<afxwin.h>

using namespace std;



#define Pi2Angle M_PI/180 // PI/180
#define MAXN 1000
#define MAXTHeight 1890
#define MAXTWidth 945
#define TEXTURE_PATH "./1.bmp"
#define SCREEN_Height 1000
#define SCREEN_Width 800
#define mLen ((SCREEN_Width + 1) * (SCREEN_Height + 1) * 3)

BYTE* DataDraw;
BYTE* DataComp;
wchar_t buffer0[50];
wchar_t buffer1[50];
wchar_t buffer2[60];
wchar_t buffer3[16];


enum CubeSide {
    LeftSide = 0,
    RightSide,
    TopSide,
    BottomSide,
    FrontSide,
    BackSide
};

struct Pixel
{
    int x;
    int y;
    COLORREF color;
    Pixel(int _x, int _y) :
        x(_x), y(_y) {};
    Pixel(int _x, int _y, COLORREF _color) :
        x(_x), y(_y), color(_color) {};
};
vector<Pixel> mPixels;

struct MRGB
{
    BYTE r;
    BYTE g;
    BYTE b;
};
//Texture[MAXN][MAXN];

class Vector {
public:
    double x, y, z, h;
    Vector(double x0 = 0, double y0 = 0, double z0 = 0, double h0 = 1) {
        x = x0;
        y = y0;
        z = z0;
        h = h0;
    }

    Vector(const Vector& v) {
        x = v.x;
        y = v.y;
        z = v.z;
        h = v.h;
    }

    ~Vector() {

    }

    Vector Reverse()const {
        Vector v;
        v.x = -x;
        v.y = -y;
        v.z = -z;
        v.h = h;

        return v;
    }

    Vector& operator= (const Vector& v) {
        if (this != &v) {
            x = v.x;
            y = v.y;
            z = v.z;
            h = v.h;
        }
        return *this;
    }

    Vector& operator+= (const Vector& v) {
        double cof = h / v.h;
        x += v.x * cof;
        y += v.y * cof;
        z += v.z * cof;

        return *this;
    }

    Vector operator* (const double f) {
        return Vector(x * f, y * f, z * f, h);
    }

    double dot(const Vector& v) {
        Vector vec1 = this->Vectorh1();
        Vector vec2 = v.Vectorh1();

        return vec1.x * vec2.x +
            vec1.y * vec2.y +
            vec1.z * vec2.z;
    }

    Vector MutiMatrix() {
        return Vector();// 需要矩阵部分
    }

    double length() const {
        Vector vec = this->Vectorh1();
        return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
    }

    Vector& Normalize() {
        double len = length();
        if (len > 0) {
            double cof = 1.0 / len;
            if (h == 0) {
                x *= cof;
                y *= cof;
                z *= cof;
            }
            else {
                x = x * cof / h;
                y = y * cof / h;
                z = z * cof / h;
            }
        }

        return *this;
    }

    Vector Vectorh1()const {
        if (h == 0) {
            return *this;
        }

        return Vector(x / h,
            y / h,
            z / h,
            1);
    }


};

Vector operator+ (const Vector& v1, const Vector& v2) {
    double cof = v1.h / v2.h;
    Vector temp;
    temp.x = v1.x + v2.x * cof;
    temp.y = v1.y + v2.y * cof;
    temp.z = v1.z + v2.z * cof;
    temp.h = v1.h;

    return temp;
}

Vector operator- (const Vector& v1, const Vector& v2) {
    double cof = v1.h / v2.h;
    Vector temp;
    temp.x = v1.x - v2.x * cof;
    temp.y = v1.y - v2.y * cof;
    temp.z = v1.z - v2.z * cof;
    temp.h = v1.h;

    return temp;
}

Vector CrossProduct(const Vector& vec1, const Vector& vec2) {
    Vector v1 = vec1.Vectorh1();
    Vector v2 = vec2.Vectorh1();

    return Vector(v1.y * v2.z - v2.y * v1.z,
        v2.x * v1.z - v1.x * v2.z,
        v1.x * v2.y - v2.x * v1.y);
}

// Line: p1p2
// Side: au+bv+cn+d = 0
Vector p_LineSide(
    const Vector& p1,
    const Vector& p2,
    double cof_u,
    double cof_v,
    double cof_n,
    double cof_costant
)
{
    double delta_u = p2.x - p1.x;
    double delta_v = p2.y - p1.y;
    double delta_n = p2.z - p1.z;
    double t = -(cof_u * p1.x + cof_v * p1.y + cof_n * p1.z + cof_costant) /
        (cof_u * delta_u + cof_v * delta_v + cof_n * delta_n);

    Vector CrossPoint;
    CrossPoint.x = p1.x + t * delta_u;
    CrossPoint.y = p1.y + t * delta_v;
    CrossPoint.z = p1.z + t * delta_n;

    return CrossPoint;
}

class TransformMatrix {
public:
    double matrix[16];

    TransformMatrix() {
        for (int i = 0; i < 16; i++) {
            matrix[i] = 0.0f;
        }
        matrix[0] = 1;
        matrix[5] = 1;
        matrix[10] = 1;
        matrix[15] = 1;
    }

    TransformMatrix(double _f11, double _f12, double _f13, double _f14,
        double _f21, double _f22, double _f23, double _f24,
        double _f31, double _f32, double _f33, double _f34,
        double _f41, double _f42, double _f43, double _f44) {
        matrix[0] = _f11;
        matrix[1] = _f12;
        matrix[2] = _f13;
        matrix[3] = _f14;
        matrix[4] = _f21;
        matrix[5] = _f22;
        matrix[6] = _f23;
        matrix[7] = _f24;
        matrix[8] = _f31;
        matrix[9] = _f32;
        matrix[10] = _f33;
        matrix[11] = _f34;
        matrix[12] = _f41;
        matrix[13] = _f42;
        matrix[14] = _f43;
        matrix[15] = _f44;
    }

    TransformMatrix(const TransformMatrix& m) {
        *this = m;
    }


    ~TransformMatrix() {

    }

    void Zero() {
        for (int i = 0; i < 16; i++) {
            matrix[i] = 0.0f;
        }
        return;
    }

    TransformMatrix& operator= (const TransformMatrix& m) {
        if (this != &m) {
            memmove(matrix, m.matrix, 16 * sizeof(double));
        }
        return *this;
    }

    TransformMatrix& operator+ (const TransformMatrix& m) {
        for (int i = 0; i < 16; i++) {
            matrix[i] += m.matrix[i];
        }
        return *this;
    }

    TransformMatrix& operator* (const TransformMatrix& m) {
        TransformMatrix M;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                M.matrix[4 * i + j] = matrix[4 * i] * m.matrix[j]
                    + matrix[4 * i + 1] * m.matrix[j + 4]
                    + matrix[4 * i + 2] * m.matrix[j + 8]
                    + matrix[4 * i + 3] * m.matrix[j + 12];
            }
        }
        return M;
    }

    Vector operator* (const Vector& v) {
        Vector vec;
        vec.x = matrix[0] * v.x
            + matrix[1] * v.y
            + matrix[2] * v.z
            + matrix[3] * v.h;
        vec.y = matrix[4] * v.x
            + matrix[5] * v.y
            + matrix[6] * v.z
            + matrix[7] * v.h;
        vec.z = matrix[8] * v.x
            + matrix[9] * v.y
            + matrix[10] * v.z
            + matrix[11] * v.h;
        vec.h = matrix[12] * v.x
            + matrix[13] * v.y
            + matrix[14] * v.z
            + matrix[15] * v.h;

        return vec;
    }

    void printTrans() {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                cout << matrix[i * 4 + j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    TransformMatrix WorldTrans(
        const Vector& position,
        const Vector& zoomVec,
        double xangle,
        double yangle,
        double zangle)
    {
        TransformMatrix M;
      
        M = M * RotateZTrans(zangle);
    
        M = M * RotateYTrans(yangle);
        M = M * RotateXTrans(xangle);
        M = M * ZoomTrans(zoomVec);
        M = M * MoveTrans(position.Reverse());
        //M.printTrans();
        return M;
    }

    TransformMatrix ObTrans(
        const Vector& VRP,
        const Vector& VPN,
        const Vector& VUP)
    {
        Vector Naxis = VPN;         //左手坐标系
        Naxis.Normalize();
        Vector Uaxis = CrossProduct(VUP, VPN);
        Uaxis.Normalize();
        Vector Vaxis = CrossProduct(Naxis, Uaxis);

        TransformMatrix M;
        M.matrix[0] = Uaxis.x;
        M.matrix[1] = Uaxis.y;
        M.matrix[2] = Uaxis.z;
        M.matrix[3] = -Uaxis.dot(VRP);
        M.matrix[4] = Vaxis.x;
        M.matrix[5] = Vaxis.y;
        M.matrix[6] = Vaxis.z;
        M.matrix[7] = -Vaxis.dot(VRP);
        M.matrix[8] = Naxis.x;
        M.matrix[9] = Naxis.y;
        M.matrix[10] = Naxis.z;
        M.matrix[11] = -Naxis.dot(VRP);

        return M;
    }

    TransformMatrix NparTrans(
        const Vector& PRP,
        double umax,
        double umin,
        double vmax,
        double vmin,
        double F,
        double B)
    {
        TransformMatrix M;
        double cof = 1 / (PRP.z - B);
        double n_min = (F - PRP.z) * cof;
        M.matrix[10] = 1 / (1 - n_min);   //统一转换为平行投影规范视见体
        M.matrix[11] = n_min / (1 - n_min);
        M.matrix[14] = -1;
        M.matrix[15] = 0;


        M = M * ZoomTrans(Vector(cof, cof, cof));
        M = M * ZoomTrans(Vector(2 * PRP.z / (umax - umin),
            2 * PRP.z / (vmax - vmin),
            1));
        M = M * SHTrans((umax + umin) / (2 * PRP.z) - PRP.x / PRP.z,
            (vmax + vmin) / (2 * PRP.z) - PRP.y / PRP.z);
        M = M * MoveTrans(PRP.Reverse());

        return M;
    }

    static TransformMatrix ProjectTrans() {
        TransformMatrix M;
        M.matrix[10] = 0;

        return M;
    }

    TransformMatrix ViewTrans(
        double height, double width, double maxz,
        double minz, double x, double y
    )
    {
        TransformMatrix M;
        M.matrix[0] = width / 2;
        M.matrix[3] = x + width / 2;
        M.matrix[5] = -height / 2;
        M.matrix[7] = y + height / 2;
        M.matrix[10] = maxz - minz;
        M.matrix[11] = minz;

        return M;
    }

    TransformMatrix MoveTrans(const Vector& moveVec) {
        TransformMatrix M;
        M.matrix[3] = moveVec.x;
        M.matrix[7] = moveVec.y;
        M.matrix[11] = moveVec.z;

        return M;
    }

    TransformMatrix ZoomTrans(const Vector& zoomVec) {
        TransformMatrix M;
        M.matrix[0] = zoomVec.x;
        M.matrix[5] = zoomVec.y;
        M.matrix[10] = zoomVec.z;

        return M;
    }

    TransformMatrix RotateXTrans(int angle) {
        TransformMatrix M;
        M.matrix[5] = cos(angle * Pi2Angle);
        M.matrix[6] = -sin(angle * Pi2Angle);
        M.matrix[9] = sin(angle * Pi2Angle);
        M.matrix[10] = cos(angle * Pi2Angle);

        return M;
    }

    TransformMatrix RotateYTrans(int angle) {
        TransformMatrix M;
        M.matrix[0] = cos(angle / 180 * M_PI);
        M.matrix[2] = sin(angle * Pi2Angle);
        M.matrix[8] = -sin(angle * Pi2Angle);
        M.matrix[10] = cos(angle / 180 * M_PI);
        return M;
    }

    TransformMatrix RotateZTrans(int angle) {
        TransformMatrix M;
        M.matrix[0] = cos(angle * Pi2Angle);
        M.matrix[1] = -sin(angle * Pi2Angle);
        M.matrix[4] = sin(angle * Pi2Angle);
        M.matrix[5] = cos(angle * Pi2Angle);

        return M;
    }

    TransformMatrix SHTrans(double cof1, double cof2) {
        TransformMatrix M;
        M.matrix[2] = cof1;
        M.matrix[6] = cof2;

        return M;
    }

    TransformMatrix RotateAxisTrans(Vector& _axis, int _angle) {
        TransformMatrix M;
        _axis.Normalize();
        double u = _axis.x;
        double v = _axis.y;
        double w = _axis.z;

        M.matrix[0] = cosf(_angle * Pi2Angle) + (u * u) * (1 - cosf(_angle * Pi2Angle));
        M.matrix[4] = u * v * (1 - cosf(_angle * Pi2Angle)) + w * sinf(_angle * Pi2Angle);
        M.matrix[8] = u * w * (1 - cosf(_angle * Pi2Angle)) - v * sinf(_angle * Pi2Angle);
        M.matrix[12] = 0;

        M.matrix[1] = u * v * (1 - cosf(_angle * Pi2Angle)) - w * sinf(_angle * Pi2Angle);
        M.matrix[5] = cosf(_angle * Pi2Angle) + v * v * (1 - cosf(_angle * Pi2Angle));
        M.matrix[9] = w * v * (1 - cosf(_angle * Pi2Angle)) + u * sinf(_angle * Pi2Angle);
        M.matrix[13] = 0;

        M.matrix[2] = u * w * (1 - cosf(_angle * Pi2Angle)) + v * sinf(_angle * Pi2Angle);
        M.matrix[6] = v * w * (1 - cosf(_angle * Pi2Angle)) - u * sinf(_angle * Pi2Angle);
        M.matrix[10] = cosf(_angle * Pi2Angle) + w * w * (1 - cosf(_angle * Pi2Angle));
        M.matrix[14] = 0;

        M.matrix[3] = 0;
        M.matrix[7] = 0;
        M.matrix[11] = 0;
        M.matrix[15] = 1;

        return M;
    }
};

class Material
{
public:
    double Ka;
    double Kd;
    double Ks;
    int n;

    Material(
        double ka = 0.0,
        double kd = 0.0,
        double ks = 0.0,
        double n_in = 0)
    {
        Ka = ka;
        Kd = kd;
        Ks = ks;
        n = n_in;
    }

    ~Material() {
    }
};

//纹理坐标
class TexCoord {
public:
    double u, v;
    TexCoord(double _u = 0.0, double _v = 0.0) {
        u = _u;
        v = _v;
    }

    TexCoord(const TexCoord& tc) {
        u = tc.u;
        v = tc.v;
    }

    TexCoord& operator=(const TexCoord& tc) {
        if (this != &tc) {
            u = tc.u;
            v = tc.v;
        }
        return *this;
    }

    ~TexCoord() {

    }
};

class Vertex {
public:
    Vector pos;         //位置
    Vector normal;      //法向量
    Material material;  //光照属性
    COLORREF color;     //颜色信息
    TexCoord tc;        //纹理坐标
    BYTE code;          //裁剪区域码
    double zdepth;      //深度值

    Vertex() {
        pos = Vector(0, 0, 0, 1);
        zdepth = 0.0;
    }

    Vertex(
        const Vector& _pos,
        const Vector& _normal,
        Material& _material,
        TexCoord& _tc)
        :pos(_pos), normal(_normal),
        material(_material), color(0),
        tc(_tc), code(0), zdepth(0.0)
    {
    }

    Vertex(const Vector& vec) {
        pos = vec;
    }

    Vertex(const Vertex& ver) {
        pos = ver.pos;
        normal = ver.normal;
        material = ver.material;
        color = ver.color;
        tc = ver.tc;
        code = ver.code;
        zdepth = ver.zdepth;
    }

    ~Vertex() {}

    Vertex& operator=(const Vertex& ver) {
        if (this != &ver) {
            pos = ver.pos;
            normal = ver.normal;
            material = ver.material;
            color = ver.color;
            tc = ver.tc;
            code = ver.code;
            zdepth = ver.zdepth;
        }

        return *this;
    }

    void Move(const Vector& vec) {
        TransformMatrix MT;
        pos = MT.MoveTrans(vec) * pos;
    }

    void Zoom(const Vector& vec) {
        TransformMatrix ZT;
        pos = ZT.ZoomTrans(vec) * pos;
    }

    void RotateX(double angle) {
        TransformMatrix R;
        pos = R.RotateXTrans(angle) * pos;
    }

    void RotateY(double angle) {
        TransformMatrix R;
        pos = R.RotateYTrans(angle) * pos;
    }

    void RotateZ(double angle) {
        TransformMatrix R;
        pos = R.RotateZTrans(angle) * pos;
    }

    void Pos_h1() {
        pos = pos.Vectorh1();
    }


};

ostream& operator<< (ostream& os, const Vertex& ver) {
    os << "坐标" << ver.pos.x << " " << ver.pos.y << " " << ver.pos.z << endl;
    os << "code" << (int)(ver.code) << endl << "color R" << (int)GetRValue(ver.color) <<" G" << (int)GetGValue(ver.color) << " B" << (int)GetBValue(ver.color) << endl;
    return os;
}
class Camera {
public:
    Vector VPN;
    Vector VRP;
    Vector VUP;
    Vector PRP;
    double umin, umax, vmin, vmax;
    double F, B;
    TransformMatrix mObTrans;
    TransformMatrix mNparTrans;

    Camera(
        Vector _VPN, Vector _VRP, Vector _VUP, Vector _PRP,
        double _umin, double _umax, double _vmin, double _vmax,
        double _F, double _B
    ) :VPN(_VPN), VRP(_VRP), VUP(_VUP), PRP(_PRP), umin(_umin),
        umax(_umax), vmin(_vmin), vmax(_vmax), F(_F), B(_B)
    {
        mObTrans = mObTrans.ObTrans(VRP, VPN, VUP);
        mNparTrans = mNparTrans.NparTrans(PRP, umax, umin, vmax, vmin, F, B);
    }

    Camera() {}

    ~Camera() {}

    Camera(const Camera& cam) {
        VPN = cam.VPN;
        VRP = cam.VRP;
        VUP = cam.VUP;
        PRP = cam.PRP;
        umin = cam.umin;
        umax = cam.umax;
        vmin = cam.vmin;
        vmax = cam.vmax;
        F = cam.F;
        B = cam.B;
        mObTrans = cam.mObTrans;
        mNparTrans = cam.mNparTrans;
    }

    Camera& operator= (const Camera& cam) {
        if (this != &cam) {
            VPN = cam.VPN;
            VRP = cam.VRP;
            VUP = cam.VUP;
            PRP = cam.PRP;
            umin = cam.umin;
            umax = cam.umax;
            vmin = cam.vmin;
            vmax = cam.vmax;
            F = cam.F;
            B = cam.B;
            mObTrans = cam.mObTrans;
            mNparTrans = cam.mNparTrans;
        }

        return *this;
    }

};

class AmbLight {
public:
    int r;
    int g;
    int b;

    AmbLight(int _r = 0, int _g = 0, int _b = 0) {
        r = _r;
        g = _g;
        b = _b;
    }

    AmbLight& operator=(const AmbLight& AL) {
        if (this != &AL) {
            r = AL.r;
            g = AL.g;
            b = AL.b;
        }

        return *this;
    }
};

class DirLight {
public:
    Vector Direction;
    int r;
    int g;
    int b;

    DirLight(const Vector& dir = Vector(0, 0, 0, 0),
        int _r = 0, int _g = 0, int _b = 0)
    {
        Direction = dir;
        r = _r;
        g = _g;
        b = _b;
    }

    ~DirLight() {}

    DirLight& operator=(const DirLight& dir) {
        if (this != &dir) {
            Direction = dir.Direction;
            r = dir.r;
            g = dir.g;
            b = dir.b;
        }

        return *this;
    }
};



class Scene {
public:
    Camera mCamera;
    // vector<Model> Models;
    AmbLight mAmbLight;
    DirLight mDirLight;

    int THeight;
    int TWidth;
    MRGB Texture[MAXTHeight][MAXTWidth];

    static Scene* mInstance;
    double x;
    double y;
    double height;
    double width;
    double minZ;
    double maxZ;
    int fps;
    bool CorrectFlag;

    Scene() {}

    ~Scene() {}

    void SetViewParam(
        double _x, double _y, double _height, double _width,
        double _minZ, double _maxZ
    )
    {
        x = _x;
        y = _y;
        height = _height;
        width = _width;
        minZ = _minZ;
        maxZ = _maxZ;
        fps = 0;
        CorrectFlag = true;
    }

    void SetCamera(const Camera& _camera) {
        mCamera = _camera;
    }

    // void AddModel(const Model& _model) {
     //    Models.push_back(_model);
   //  }

    static Scene* GetInstance() {
        if (mInstance == NULL) {
            mInstance = new Scene();
        }
        return mInstance;
    }

    void SetAmbLight(const AmbLight& _ambt) {
        mAmbLight = _ambt;
    }

    void SetDirLight(const DirLight& _paral) {
        mDirLight = _paral;
    }

    void ReadTexture() {
        BITMAPFILEHEADER fileHeader;
        BITMAPINFOHEADER infoHeader;
        FILE* pfin = fopen(TEXTURE_PATH, "rb");
        //Read Bitmap file header
        fread(&fileHeader, sizeof(BITMAPFILEHEADER), 1, pfin);
        //Read Bitmap info header
        fread(&infoHeader, sizeof(BITMAPINFOHEADER), 1, pfin);
        //24 Bit color
        THeight = abs(infoHeader.biHeight);
        TWidth = abs(infoHeader.biWidth);
        int size = THeight * TWidth;
        fread(Texture, sizeof(MRGB), size, pfin);
        //cout << THeight << endl << TWidth << endl;
        /*for (int i = 0; i < THeight; i++) {
            for (int j = 0; j < TWidth; j++) {
                if ((int)Texture[i][j].r == 205)
                    cout << (int)Texture[i][j].r << " " << (int)Texture[i][j].r << " " << (int)Texture[i][j].r << "   ";
            }
            cout << endl;
        }
        */
    }

    //void Draw() {
    //    for (auto iter = Models.begin(); iter != Models.end(); iter++) {
    //        iter->Draw();
    //    }
   // }

};
Scene* Scene::mInstance = NULL;


class sTriangle {
public:
    Vertex ver[3];
    sTriangle() {}
    sTriangle(
        const Vertex& ver0,
        const Vertex& ver1,
        const Vertex& ver2)
    {
        ver[0] = ver0;
        ver[1] = ver1;
        ver[2] = ver2;
    }
    ~sTriangle() {}

    void printSTri()const {
        for (int i = 0; i < 3; i++) {
            cout << ver[i] << endl;
        }
    }
};

class Triangle {
public:
    Vertex pos[3];


    Triangle() {
        memset(pos, 0, sizeof(pos));
        Init();
    }

    Triangle(
        const Vertex& pos0,
        const Vertex& pos1,
        const Vertex& pos2)
    {
        pos[0] = pos0;
        pos[1] = pos1;
        pos[2] = pos2;

        Init();
    }

    ~Triangle() {}

    void Init() {

    }

    Triangle(const Triangle& tri) {
        *this = tri;
        Init();
    }

    Triangle& operator=(const Triangle& tri) {
        if (this != &tri) {
            pos[0] = tri.pos[0];
            pos[1] = tri.pos[1];
            pos[2] = tri.pos[2];
        }

        Init();
        return *this;
    }

    COLORREF LightModel(const Vertex& ver) {
        AmbLight Ia = Scene::GetInstance()->mAmbLight;
        DirLight Ip = Scene::GetInstance()->mDirLight;

        Vector L = Ip.Direction.Reverse();
        L.Normalize();
        Vector N = ver.normal;
        N.Normalize();

        Camera cam = Scene::GetInstance()->mCamera;
        Vector V = cam.VRP - ver.pos;
        V.Normalize();
        Vector H = (L + V) * (1 / 2);
        H.Normalize();

        int Ir;
        int Ig;
        int Ib;

        Material m = ver.material;
        //漫反射及镜面反射因子
        double factor = m.Kd * (L.dot(N)) + m.Ks * pow((H.dot(N)), m.n);
        factor = factor < 0 ? 0 : factor;

        Ir = m.Ka * Ia.r + Ip.r * factor;
        if (Ir < 0) Ir = 0;
        else if (Ir > 255) Ir = 255;
        Ig = m.Ka * Ia.g + Ip.g * factor;
        if (Ig < 0) Ig = 0;
        else if (Ig > 255) Ig = 255;
        Ib = m.Ka * Ia.b + Ip.b * factor;
        if (Ib < 0) Ib = 0;
        else if (Ib > 255) Ib = 255;

        return COLORREF(RGB(Ir, Ig, Ib));
    }

    void SwapVertex(Vertex& ver1, Vertex& ver2) {
        Vertex temp = ver2;
        ver2 = ver1;
        ver1 = temp;
        return;
    }

    COLORREF ColorCrossX(const Vertex& p0, const Vertex& p1, double x) {
        double newR = GetRValue(p0.color) +
            (GetRValue(p1.color) - GetRValue(p0.color)) /
            (p1.pos.x - p0.pos.x) * (x - p0.pos.x);
        double newG = GetGValue(p0.color) +
            (GetGValue(p1.color) - GetGValue(p0.color)) /
            (p1.pos.x - p0.pos.x) * (x - p0.pos.x);
        double newB = GetBValue(p0.color) +
            (GetBValue(p1.color) - GetBValue(p0.color)) /
            (p1.pos.x - p0.pos.x) * (x - p0.pos.x);
        return RGB(newR, newG, newB);
    }

    COLORREF ColorCrossY(const Vertex& p0, const Vertex& p1, double y) {
        double newR = GetRValue(p0.color) +
            (GetRValue(p1.color) - GetRValue(p0.color)) /
            (p1.pos.y - p0.pos.y) * (y - p0.pos.y);
        double newG = GetGValue(p0.color) +
            (GetGValue(p1.color) - GetGValue(p0.color)) /
            (p1.pos.y - p0.pos.y) * (y - p0.pos.y);
        double newB = GetBValue(p0.color) +
            (GetBValue(p1.color) - GetBValue(p0.color)) /
            (p1.pos.y - p0.pos.y) * (y - p0.pos.y);
        return RGB(newR, newG, newB);
    }

    COLORREF ColorCrossZ(const Vertex& p0, const Vertex& p1, double z) {
        double newR = GetRValue(p0.color) +
            (GetRValue(p1.color) - GetRValue(p0.color)) /
            (p1.pos.z - p0.pos.z) * (z - p0.pos.z);
        double newG = GetGValue(p0.color) +
            (GetGValue(p1.color) - GetGValue(p0.color)) /
            (p1.pos.z - p0.pos.z) * (z - p0.pos.z);
        double newB = GetBValue(p0.color) +
            (GetBValue(p1.color) - GetBValue(p0.color)) /
            (p1.pos.z - p0.pos.z) * (z - p0.pos.z);
        return RGB(newR, newG, newB);
    }

    //输入一个三角形和一个面，三角形可能被裁剪成0-2个三角形
    //此时保证三角形和平面一定相交
    void Cut(
        const sTriangle& _tri,
        CubeSide _side,
        vector<sTriangle>& _out
    )
    {
        switch (_side)
        {
            // v=1
        case TopSide:
        {
            //将三个顶点按上下顺序排序
            Vertex ver0 = _tri.ver[0];
            Vertex ver1 = _tri.ver[1];
            Vertex ver2 = _tri.ver[2];

            if (ver0.pos.y < ver2.pos.y)
                SwapVertex(ver0, ver2);
            if (ver0.pos.y < ver1.pos.y)
                SwapVertex(ver0, ver1);
            if (ver1.pos.y < ver2.pos.y)
                SwapVertex(ver1, ver2);

            //三点一线则无法构成三角形
            if (abs(ver0.pos.y - ver2.pos.y) < 0.0000001f) {
                return;
            }

            //此时三角形与v=1相交情况为： 
            //  ver0ver2必定相交，ver0ver1还是ver1ver2与v=1相交需要分类讨论
            Vertex v0;
            v0.pos = p_LineSide(ver0.pos, ver2.pos, 0, 1, 0, -1);
            v0.normal = ver0.normal;
            //该点的颜色做线性插值
            v0.color = ColorCrossY(ver0, ver2, 1);

            //v0纹理坐标做透视修正
            // 还没看*****************
            double s = (1.0 - ver0.pos.y) /
                (ver2.pos.y - ver0.pos.y);
            double Z1 = 0.0;
            double Z2 = 0.0;
            double Zt = 0.0;
            double t = s;
            Z1 = ver0.zdepth;
            Z2 = ver2.zdepth;
            if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
            if (Zt != 0) Zt = 1 / Zt;
            if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
            v0.zdepth = Zt;
            v0.tc.u = ver0.tc.u +
                t * (ver2.tc.u - ver0.tc.u);
            v0.tc.v = ver0.tc.v +
                t * (ver2.tc.v - ver0.tc.v);

            //ver0在平面上 ver1 ver2在平面下 即ver0ver1与平面相交
            if (ver0.pos.y > 1.0f && ver1.pos.y < 1.0f) {
                Vertex v1;
                v1.pos = p_LineSide(ver0.pos, ver1.pos, 0, 1, 0, -1);
                v1.normal = ver0.normal;
                v1.color = ColorCrossY(ver0, ver1, 1);

                //v1纹理坐标做透视修正
                // 还没看*****************
                double s = (1.0 - ver0.pos.y) /
                    (ver1.pos.y - ver0.pos.y);
                double Z1 = 0.0;
                double Z2 = 0.0;
                double Zt = 0.0;
                double t = s;
                Z1 = ver0.zdepth;
                Z2 = ver1.zdepth;
                if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
                if (Zt != 0) Zt = 1 / Zt;
                if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
                v1.zdepth = Zt;
                v1.tc.u = ver0.tc.u +
                    t * (ver1.tc.u - ver0.tc.u);
                v1.tc.v = ver0.tc.v +
                    t * (ver1.tc.v - ver0.tc.v);

                _out.push_back(sTriangle(v0, v1, ver2));
                _out.push_back(sTriangle(ver1, ver2, v1));
            }
            //ver0 ver1 在平面上 ver2在平面下 即ver1ver2与平面相交
            else if (ver0.pos.y > 1.0f && ver1.pos.y >= 1.0f) {
                Vertex v2;
                v2.pos = p_LineSide(ver1.pos, ver2.pos, 0, 1, 0, -1);
                v2.normal = ver1.normal;
                v2.color = ColorCrossY(ver1, ver2, 1);

                //v2纹理坐标做透视修正
                // 还没看*****************
                double s = (1.0 - ver1.pos.y) /
                    (ver2.pos.y - ver1.pos.y);
                double Z1 = 0.0;
                double Z2 = 0.0;
                double Zt = 0.0;
                double t = s;
                Z1 = ver1.zdepth;
                Z2 = ver2.zdepth;
                if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
                if (Zt != 0) Zt = 1 / Zt;
                if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
                v2.zdepth = Zt;
                v2.tc.u = ver1.tc.u +
                    t * (ver2.tc.u - ver1.tc.u);
                v2.tc.v = ver1.tc.v +
                    t * (ver2.tc.v - ver1.tc.v);

                _out.push_back(sTriangle(v0, v2, ver2));
            }
            break;
        }
        //v=-1
        case BottomSide:
        {
            //将三个顶点按上下顺序排序
            Vertex ver0 = _tri.ver[0];
            Vertex ver1 = _tri.ver[1];
            Vertex ver2 = _tri.ver[2];

            if (ver0.pos.y < ver2.pos.y)
                SwapVertex(ver0, ver2);
            if (ver0.pos.y < ver1.pos.y)
                SwapVertex(ver0, ver1);
            if (ver1.pos.y < ver2.pos.y)
                SwapVertex(ver1, ver2);

            //三点一线则无法构成三角形
            if (abs(ver0.pos.y - ver2.pos.y) < 0.0000001f) {
                return;
            }

            //此时三角形与v=-1相交情况为： 
            //  ver0ver2必定相交，ver0ver1还是ver1ver2与v=-1相交需要分类讨论
            Vertex v0;
            v0.pos = p_LineSide(ver0.pos, ver2.pos, 0, 1, 0, 1);
            v0.normal = ver0.normal;
            //该点的颜色做线性插值
            v0.color = ColorCrossY(ver0, ver2, -1);

            //v0纹理坐标做透视修正
            // 还没看*****************
            double s = (-1.0 - ver0.pos.y) /
                (ver2.pos.y - ver0.pos.y);
            double Z1 = 0.0;
            double Z2 = 0.0;
            double Zt = 0.0;
            double t = s;
            Z1 = ver0.zdepth;
            Z2 = ver2.zdepth;
            if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
            if (Zt != 0) Zt = 1 / Zt;
            if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
            v0.zdepth = Zt;
            v0.tc.u = ver0.tc.u +
                t * (ver2.tc.u - ver0.tc.u);
            v0.tc.v = ver0.tc.v +
                t * (ver2.tc.v - ver0.tc.v);


            //ver0在平面上 ver1 ver2在平面下 即ver0ver1与平面相交
            if (ver0.pos.y > -1.0f && ver1.pos.y < -1.0f) {
                Vertex v1;
                v1.pos = p_LineSide(ver0.pos, ver1.pos, 0, 1, 0, 1);
                v1.normal = ver0.normal;
                v1.color = ColorCrossY(ver0, ver1, -1);

                //v1纹理坐标做透视修正
                // 还没看*****************
                double s = (-1.0 - ver0.pos.y) /
                    (ver1.pos.y - ver0.pos.y);
                double Z1 = 0.0;
                double Z2 = 0.0;
                double Zt = 0.0;
                double t = s;
                Z1 = ver0.zdepth;
                Z2 = ver1.zdepth;
                if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
                if (Zt != 0) Zt = 1 / Zt;
                if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
                v1.zdepth = Zt;
                v1.tc.u = ver0.tc.u +
                    t * (ver1.tc.u - ver0.tc.u);
                v1.tc.v = ver0.tc.v +
                    t * (ver1.tc.v - ver0.tc.v);

                _out.push_back(sTriangle(v0, v1, ver0));
            }
            //ver0 ver1 在平面上 ver2在平面下 即ver1ver2与平面相交
            else if (ver0.pos.y > -1.0f && ver1.pos.y >= -1.0f) {
                Vertex v2;
                v2.pos = p_LineSide(ver1.pos, ver2.pos, 0, 1, 0, 1);
                v2.normal = ver1.normal;
                v2.color = ColorCrossY(ver1, ver2, -1);

                //v2纹理坐标做透视修正
                // 还没看*****************
                double s = (-1.0 - ver1.pos.y) /
                    (ver2.pos.y - ver1.pos.y);
                double Z1 = 0.0;
                double Z2 = 0.0;
                double Zt = 0.0;
                double t = s;
                Z1 = ver1.zdepth;
                Z2 = ver2.zdepth;
                if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
                if (Zt != 0) Zt = 1 / Zt;
                if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
                v2.zdepth = Zt;
                v2.tc.u = ver1.tc.u +
                    t * (ver2.tc.u - ver1.tc.u);
                v2.tc.v = ver1.tc.v +
                    t * (ver2.tc.v - ver1.tc.v);

                _out.push_back(sTriangle(v0, v2, ver0));
                _out.push_back(sTriangle(ver0, ver1, v2));
            }
            break;
        }
        //x = -1
        case LeftSide:
        {
            //将顶点按x排序
            Vertex ver0 = _tri.ver[0];
            Vertex ver1 = _tri.ver[1];
            Vertex ver2 = _tri.ver[2];

            if (ver0.pos.x > ver2.pos.x)
                SwapVertex(ver0, ver2);
            if (ver0.pos.x > ver1.pos.x)
                SwapVertex(ver0, ver1);
            if (ver1.pos.x > ver2.pos.x)
                SwapVertex(ver1, ver2);

            if (abs(ver0.pos.x - ver2.pos.x) < 0.0000001f)
                return;

            Vertex v0;
            v0.pos = p_LineSide(ver0.pos, ver2.pos, 1, 0, 0, 1);
            v0.normal = ver0.normal;
            v0.color = ColorCrossX(ver0, ver2, -1);

            //透视矫正
            double s = (-1.0 - ver0.pos.x) /
                (ver2.pos.x - ver0.pos.x);
            double Z1 = 0.0;
            double Z2 = 0.0;
            double Zt = 0.0;
            double t = s;
            Z1 = ver0.zdepth;
            Z2 = ver2.zdepth;
            if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
            if (Zt != 0) Zt = 1 / Zt;
            if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
            v0.zdepth = Zt;
            v0.tc.u = ver0.tc.u +
                t * (ver2.tc.u - ver0.tc.u);
            v0.tc.v = ver0.tc.v +
                t * (ver2.tc.v - ver0.tc.v);


            //ver0ver1与x=-1相交
            if (ver0.pos.x < -1.0f && ver1.pos.x > -1.0f) {
                Vertex v1;
                v1.pos = p_LineSide(ver0.pos, ver1.pos, 1, 0, 0, 1);
                v1.normal = ver0.normal;
                v1.color = ColorCrossX(ver0, ver1, -1);

                //透视矫正
                double s = (-1.0 - ver0.pos.x) /
                    (ver1.pos.x - ver0.pos.x);
                double Z1 = 0.0;
                double Z2 = 0.0;
                double Zt = 0.0;
                double t = s;
                Z1 = ver0.zdepth;
                Z2 = ver1.zdepth;
                if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
                if (Zt != 0) Zt = 1 / Zt;
                if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
                v1.zdepth = Zt;
                v1.tc.u = ver0.tc.u +
                    t * (ver1.tc.u - ver0.tc.u);
                v1.tc.v = ver0.tc.v +
                    t * (ver1.tc.v - ver0.tc.v);

                _out.push_back(sTriangle(v0, v1, ver2));
                _out.push_back(sTriangle(ver1, ver2, v1));
            }
            else if (ver0.pos.x < -1.0f && ver1.pos.x <= -1.0f) {
                Vertex v2;
                v2.pos = p_LineSide(ver1.pos, ver2.pos, 1, 0, 0, 1);
                v2.normal = ver1.normal;
                v2.color = ColorCrossX(ver1, ver2, -1);

                //透视矫正
                double s = (-1.0 - ver1.pos.x) /
                    (ver2.pos.x - ver1.pos.x);
                double Z1 = 0.0;
                double Z2 = 0.0;
                double Zt = 0.0;
                double t = s;
                Z1 = ver1.zdepth;
                Z2 = ver2.zdepth;
                if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
                if (Zt != 0) Zt = 1 / Zt;
                if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
                v2.zdepth = Zt;
                v2.tc.u = ver1.tc.u +
                    t * (ver2.tc.u - ver1.tc.u);
                v2.tc.v = ver1.tc.v +
                    t * (ver2.tc.v - ver1.tc.v);

                _out.push_back(sTriangle(v0, v2, ver2));
            }
            break;
        }
        // x=1
        case RightSide:
        {
            //将顶点按x排序
            Vertex ver0 = _tri.ver[0];
            Vertex ver1 = _tri.ver[1];
            Vertex ver2 = _tri.ver[2];

            if (ver0.pos.x > ver2.pos.x)
                SwapVertex(ver0, ver2);
            if (ver0.pos.x > ver1.pos.x)
                SwapVertex(ver0, ver1);
            if (ver1.pos.x > ver2.pos.x)
                SwapVertex(ver1, ver2);

            if (abs(ver0.pos.x - ver2.pos.x) < 0.0000001f)
                return;

            Vertex v0;
            v0.pos = p_LineSide(ver0.pos, ver2.pos, 1, 0, 0, -1);
            v0.normal = ver0.normal;
            v0.color = ColorCrossX(ver0, ver2, 1);

            //透视矫正
            double s = (1.0 - ver0.pos.x) /
                (ver2.pos.x - ver0.pos.x);
            double Z1 = 0.0;
            double Z2 = 0.0;
            double Zt = 0.0;
            double t = s;
            Z1 = ver0.zdepth;
            Z2 = ver2.zdepth;
            if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
            if (Zt != 0) Zt = 1 / Zt;
            if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
            v0.zdepth = Zt;
            v0.tc.u = ver0.tc.u +
                t * (ver2.tc.u - ver0.tc.u);
            v0.tc.v = ver0.tc.v +
                t * (ver2.tc.v - ver0.tc.v);

            //ver0ver1与x=1相交
            if (ver0.pos.x < 1.0f && ver1.pos.x > 1.0f) {
                Vertex v1;
                v1.pos = p_LineSide(ver0.pos, ver1.pos, 1, 0, 0, -1);
                v1.normal = ver0.normal;
                v1.color = ColorCrossX(ver0, ver1, 1);

                //透视矫正
                double s = (1.0 - ver0.pos.x) /
                    (ver1.pos.x - ver0.pos.x);
                double Z1 = 0.0;
                double Z2 = 0.0;
                double Zt = 0.0;
                double t = s;
                Z1 = ver0.zdepth;
                Z2 = ver1.zdepth;
                if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
                if (Zt != 0) Zt = 1 / Zt;
                if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
                v1.zdepth = Zt;
                v1.tc.u = ver0.tc.u +
                    t * (ver1.tc.u - ver0.tc.u);
                v1.tc.v = ver0.tc.v +
                    t * (ver1.tc.v - ver0.tc.v);

                _out.push_back(sTriangle(v0, v1, ver0));
            }
            else if (ver0.pos.x < 1.0f && ver1.pos.x <= 1.0f) {
                Vertex v2;
                v2.pos = p_LineSide(ver1.pos, ver2.pos, 1, 0, 0, -1);
                v2.normal = ver1.normal;
                v2.color = ColorCrossX(ver1, ver2, 1);

                //透视矫正
                double s = (1.0 - ver1.pos.x) /
                    (ver2.pos.x - ver1.pos.x);
                double Z1 = 0.0;
                double Z2 = 0.0;
                double Zt = 0.0;
                double t = s;
                Z1 = ver1.zdepth;
                Z2 = ver2.zdepth;
                if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
                if (Zt != 0) Zt = 1 / Zt;
                if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
                v2.zdepth = Zt;
                v2.tc.u = ver1.tc.u +
                    t * (ver2.tc.u - ver1.tc.u);
                v2.tc.v = ver1.tc.v +
                    t * (ver2.tc.v - ver1.tc.v);

                _out.push_back(sTriangle(ver0, ver1, v2));
                _out.push_back(sTriangle(v0, v2, ver0));
            }
            break;
        }
        //z=-1
        case FrontSide:
        {
            Vertex ver0 = _tri.ver[0];
            Vertex ver1 = _tri.ver[1];
            Vertex ver2 = _tri.ver[2];

            //按z由由小到大排序
            if (ver0.pos.z > ver2.pos.z)
                SwapVertex(ver0, ver2);
            if (ver0.pos.z > ver1.pos.z)
                SwapVertex(ver0, ver1);
            if (ver1.pos.z > ver2.pos.z)
                SwapVertex(ver1, ver2);

            if (abs(ver0.pos.z - ver2.pos.z) < 0.00000001f)
                return;
            cout << ver0 << endl << ver1 << endl << ver2 << endl;
            Vertex v0;
            v0.pos = p_LineSide(ver0.pos, ver2.pos, 0, 0, 1, 1);
            v0.normal = ver0.normal;
            v0.color = ColorCrossZ(ver0, ver2, -1);

            cout << v0 << endl;
            //透视矫正
            double s = (-1.0 - ver0.pos.z) / (ver2.pos.z - ver0.pos.z);
            double Z1 = 0.0;
            double Z2 = 0.0;
            double Zt = 0.0;
            double t = s;
            Z1 = ver0.zdepth;
            Z2 = ver2.zdepth;
            if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
            if (Zt != 0) Zt = 1 / Zt;
            if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
            v0.zdepth = Zt;
            v0.tc.u = ver0.tc.u +
                t * (ver2.tc.u - ver0.tc.u);
            v0.tc.v = ver0.tc.v +
                t * (ver2.tc.v - ver0.tc.v);

            //ver0ver1与z=-1相交
            if (ver0.pos.z < -1.0f && ver1.pos.z > -1.0f) {
                Vertex v1;
                v1.pos = p_LineSide(ver0.pos, ver1.pos, 0, 0, 1, 1);
                v1.normal = ver0.normal;
                v1.color = ColorCrossZ(ver0, ver1, -1);

                //透视修正
                double s = (-1.0 - ver0.pos.z) / (ver1.pos.z - ver0.pos.z);
                double Z1 = 0.0;
                double Z2 = 0.0;
                double Zt = 0.0;
                double t = s;
                Z1 = ver0.zdepth;
                Z2 = ver1.zdepth;
                if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
                if (Zt != 0) Zt = 1 / Zt;
                if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
                v1.zdepth = Zt;
                v1.tc.u = ver0.tc.u +
                    t * (ver1.tc.u - ver0.tc.u);
                v1.tc.v = ver0.tc.v +
                    t * (ver1.tc.v - ver0.tc.v);
                _out.push_back(sTriangle(ver1, ver2, v1));
                _out.push_back(sTriangle(v0, v1, ver2));
            }
            else if (ver0.pos.z < -1.0f && ver1.pos.z <= -1.0f)
            {
                Vertex v2;
                v2.pos = p_LineSide(ver1.pos, ver2.pos, 0, 0, 1, 1);
                v2.normal = ver1.normal;
                v2.color = ColorCrossZ(ver1, ver2, -1);

                //透视矫正
                double s = (-1.0 - ver1.pos.z) / (ver2.pos.z - ver1.pos.z);
                double Z1 = 0.0;
                double Z2 = 0.0;
                double Zt = 0.0;
                double t = s;
                Z1 = ver1.zdepth;
                Z2 = ver2.zdepth;
                if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
                if (Zt != 0) Zt = 1 / Zt;
                if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
                v2.zdepth = Zt;
                v2.tc.u = ver1.tc.u +
                    t * (ver2.tc.u - ver1.tc.u);
                v2.tc.v = ver1.tc.v +
                    t * (ver2.tc.v - ver1.tc.v);

                _out.push_back(sTriangle(v0, v2, ver2));
            }
            break;
        }
        //z=0
        case BackSide:
        {
            Vertex ver0 = _tri.ver[0];
            Vertex ver1 = _tri.ver[1];
            Vertex ver2 = _tri.ver[2];

            //按z由由小到大排序
            if (ver0.pos.z > ver2.pos.z)
                SwapVertex(ver0, ver2);
            if (ver0.pos.z > ver1.pos.z)
                SwapVertex(ver0, ver1);
            if (ver1.pos.z > ver2.pos.z)
                SwapVertex(ver1, ver2);

            if (abs(ver0.pos.z - ver2.pos.z) < 0.00000001f)
                return;

            Vertex v0;
            v0.pos = p_LineSide(ver0.pos, ver2.pos, 0, 0, 1, 0);
            v0.normal = ver0.normal;
            v0.color = ColorCrossZ(ver0, ver2, 0);

            //透视矫正
            double s = (0.0 - ver0.pos.z) / (ver2.pos.z - ver0.pos.z);
            double Z1 = 0.0;
            double Z2 = 0.0;
            double Zt = 0.0;
            double t = s;
            Z1 = ver0.zdepth;
            Z2 = ver2.zdepth;
            if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
            if (Zt != 0) Zt = 1 / Zt;
            if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
            v0.zdepth = Zt;
            v0.tc.u = ver0.tc.u +
                t * (ver2.tc.u - ver0.tc.u);
            v0.tc.v = ver0.tc.v +
                t * (ver2.tc.v - ver0.tc.v);

            //ver0ver1与z=0相交
            if (ver0.pos.z < 0.0000001f && ver1.pos.z > 0.0000001f) {
                Vertex v1;
                v1.pos = p_LineSide(ver0.pos, ver1.pos, 0, 0, 1, 0);
                v1.normal = ver0.normal;
                v1.color = ColorCrossZ(ver0, ver1, 0);

                //透视修正
                double s = (0.0 - ver0.pos.z) / (ver1.pos.z - ver0.pos.z);
                double Z1 = 0.0;
                double Z2 = 0.0;
                double Zt = 0.0;
                double t = s;
                Z1 = ver0.zdepth;
                Z2 = ver1.zdepth;
                if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
                if (Zt != 0) Zt = 1 / Zt;
                if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
                v1.zdepth = Zt;
                v1.tc.u = ver0.tc.u +
                    t * (ver1.tc.u - ver0.tc.u);
                v1.tc.v = ver0.tc.v +
                    t * (ver1.tc.v - ver0.tc.v);
                _out.push_back(sTriangle(v0, v1, ver0));
            }
            else if (ver0.pos.z < 0.0000001f && ver1.pos.z <= 0.0000001f)
            {
                Vertex v2;
                v2.pos = p_LineSide(ver1.pos, ver2.pos, 0, 0, 1, 0);
                v2.normal = ver1.normal;
                v2.color = ColorCrossZ(ver1, ver2, 0);

                //透视矫正
                double s = (0.0 - ver1.pos.z) / (ver2.pos.z - ver1.pos.z);
                double Z1 = 0.0;
                double Z2 = 0.0;
                double Zt = 0.0;
                double t = s;
                Z1 = ver1.zdepth;
                Z2 = ver2.zdepth;
                if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
                if (Zt != 0) Zt = 1 / Zt;
                if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
                v2.zdepth = Zt;
                v2.tc.u = ver1.tc.u +
                    t * (ver2.tc.u - ver1.tc.u);
                v2.tc.v = ver1.tc.v +
                    t * (ver2.tc.v - ver1.tc.v);
                _out.push_back(sTriangle(v0, v2, ver0));
                _out.push_back(sTriangle(ver0, ver1, v2));
            }
            break;
        }
        default:
            break;
        }
    }

    //C0C1C2C3C4C5对应6个面的内外部情况
    void CalculateCode(Vertex& ver) {
        BYTE AreaCode = 0;
        if (ver.pos.x < -1.0f)   AreaCode |= 1;
        if (ver.pos.x > 1.0f)    AreaCode |= 2;
        if (ver.pos.y < -1.0f)   AreaCode |= 4;
        if (ver.pos.y > 1.0f)    AreaCode |= 8;
        if (ver.pos.z < -1.0f)   AreaCode |= 16;
        if (ver.pos.z > 0.0000001f) AreaCode |= 32;

        ver.code = AreaCode;
        return;
    }

    //C5 C4 C3 C2 C1 C0
    //z=0/-1 y=1/-1 x=1/-1
    void SideCrossOfL(
        const Vertex& ver0,
        const Vertex& ver1,
        vector<CubeSide>& sides
    )
    {
        int side = 0; // x=-1
      //  if (((ver0.code & ver1.code) == 0) || ((ver0.code | ver1.code) != 0))
        //    return;

        BYTE c = (ver0.code ^ ver1.code);
        while (c) {
            if ((c & 1) == 1) {
                sides.push_back(CubeSide(side));
            }
            side++;
            c >>= 1;
        }
    }

    void SideCrossofTri(
        const sTriangle& tri,
        queue<CubeSide>& sides
    )
    {
        vector<CubeSide> temp;
        SideCrossOfL(tri.ver[0], tri.ver[1], temp);
        SideCrossOfL(tri.ver[0], tri.ver[2], temp);
        SideCrossOfL(tri.ver[1], tri.ver[2], temp);
 
        for (int i = 0; i < 6; i++) {
            auto iter = find(temp.begin(), temp.end(), (CubeSide)i);
            if (iter != temp.end()) {
                sides.push((CubeSide)i);
            }
        }
    }

    //规范视见体裁剪三角形
    void Cut(
        const sTriangle& tri,
        vector<sTriangle>& outTris
    )
    {

        if ((tri.ver[0].code & tri.ver[1].code & tri.ver[2].code) > 0)
            return;
        else if ((tri.ver[0].code | tri.ver[1].code | tri.ver[2].code) == 0)
            outTris.push_back(tri);
        else {
            queue<CubeSide> Sides;
            SideCrossofTri(tri, Sides);
            //cout << (CubeSide)Sides.front();
            vector<sTriangle> Tris_in;
            Tris_in.push_back(tri);

            while (!Sides.empty()) {
                auto Side = Sides.front();
                Sides.pop();

                for (auto iter = Tris_in.begin();
                    iter != Tris_in.end(); iter++)
                {
                    sTriangle CurrentTri = *iter;
                    //CurrentTri.printSTri();
                    //cout << Side << endl;
                    Cut(CurrentTri, Side, outTris);
                }
                //outTris.front().printSTri();
                if (Sides.empty());
                    break;
                Tris_in.clear();

                for (auto iter = outTris.begin();
                    iter != outTris.end(); iter++)
                {
                    Tris_in.push_back(*iter);
                }
                outTris.clear();
            }
            //cout << outTris.size() << endl;
        }
    }

    void DDALine(double x0, double x1, double y,
        double r0, double g0, double b0, double r1, double g1, double b1,
        double z0, double u0, double v0, double z1, double u1, double v1)
    {
        double dr = 0.0;
        double dg = 0.0;
        double db = 0.0;

        if (abs(x1 - x0) > 0.000001f) {
            dr = (r1 - r0) / (x1 - x0);
            dg = (g1 - g0) / (x1 - x0);
            db = (b1 - b0) / (x1 - x0);
        }

        double R = r0;
        double G = g0;
        double B = b0;

        for (double x = x0; x <= x1; x += 0.8f) {
            if (abs(x1 - x0) > 0.000001f) {
                double s = (x - x0) / (x1 - x0);
                double Zt = 0;
                double t = s;
                double Z0 = z0;
                double Z2 = z1;
                if (Z0 != 0.0 && Z2 != 0.0) Zt = 1 / Z0 + s * (1 / Z2 - 1 / Z0);
                if (Zt != 0) Zt = 1 / Zt;
                if (Z2 != Z0) t = (Zt - Z0) / (Z2 - Z0);

                double u = u0 + t * (u1 - u0);
                double v = v0 + t * (v1 - v0);
                if (v < 0.0) v = 0.0;
                else if (v > 1.0f) v = 1.0;
                if (u < 0.0) u = 0.0;
                else if (u > 1.0f) u = 1.0;

                int THeight = Scene::GetInstance()->THeight;
                int TWidth = Scene::GetInstance()->TWidth;

                MRGB textcolor = Scene::GetInstance()->Texture[(int)((THeight - 1) * u)][(int)((TWidth - 1) * v)];

                int r = R * ((double)(textcolor.r) / 255.0);
                int g = G * ((double)(textcolor.g) / 255.0);
                int b = B * ((double)(textcolor.b) / 255.0);
                //cout << x << " " << y << endl << (int)textcolor.r << " " << (int)textcolor.g << " " << (int)textcolor.b << endl;
                //putpixel(x, y, RGB(r, g, b));
                mPixels.push_back(Pixel(x, y, RGB(r, g, b)));
            }
            else {
                double u = u0;
                double v = v0;
                if (v < 0.0) v = 0.0;
                else if (v > 1.0f) v = 1.0;
                if (u < 0.0) u = 0.0;
                else if (u > 1.0f) u = 1.0;

                int THeight = Scene::GetInstance()->THeight;
                int TWidth = Scene::GetInstance()->TWidth;

                MRGB textcolor = Scene::GetInstance()->Texture[(int)((THeight - 1) * u)][(int)((TWidth - 1) * v)];
                if ((int)textcolor.r == 205) {
                    textcolor.r = 255;
                    textcolor.g = 255;
                    textcolor.b = 255;
                }
                int r = R * ((double)(textcolor.r) / 255.0);
                int g = G * ((double)(textcolor.g) / 255.0);
                int b = B * ((double)(textcolor.b) / 255.0);

                //putpixel(x, y, RGB(r, g, b));
                mPixels.push_back(Pixel(x, y, RGB(r, g, b)));
                //cout << x << " " << y << endl << textcolor.r << " " << textcolor.g << " " << textcolor.b << endl;
            }
            R += dr;
            G += dg;
            B += db;
        }
    }

    void DrawTopTriangle(Vertex& _ver0, Vertex& _ver1, Vertex& _ver2) {
        //判断输入的三角形
        if (abs(_ver0.pos.y - _ver1.pos.y) < 0.000001f)
        {/*empty*/
        }
        else if (abs(_ver0.pos.y - _ver2.pos.y) < 0.000001f)
        {
            SwapVertex(_ver2, _ver1);
        }
        else if (abs(_ver1.pos.y - _ver2.pos.y) < 0.000001f)
        {
            SwapVertex(_ver0, _ver2);
        }
        else
        { //不是平顶三角形
            return;
        }
        if (_ver1.pos.x < _ver0.pos.x)
        {
            SwapVertex(_ver1, _ver0);
        }
        else if (abs(_ver1.pos.x - _ver0.pos.x) < 0.000001f)
        { //不是三角形
            return;
        }

        double _x0 = _ver0.pos.x;
        double _x1 = _ver1.pos.x;
        double _x2 = _ver2.pos.x;
        double _y0 = _ver0.pos.y;
        double _y1 = _ver1.pos.y;
        double _y2 = _ver2.pos.y;
        COLORREF _color0 = _ver0.color;
        COLORREF _color1 = _ver1.color;
        COLORREF _color2 = _ver2.color;

        //计算左右误差
        double dxy_left = (_x2 - _x0) * 1.0 / (_y2 - _y0);
        double dxy_right = (_x1 - _x2) * 1.0 / (_y1 - _y2);
        double dr_left = (GetRValue(_color2) - GetRValue(_color0)) * 1.0 / (_y2 - _y0);
        double dr_right = (GetRValue(_color1) - GetRValue(_color2)) * 1.0 / (_y1 - _y2);
        double dg_left = (GetGValue(_color2) - GetGValue(_color0)) * 1.0 / (_y2 - _y0);
        double dg_right = (GetGValue(_color1) - GetGValue(_color2)) * 1.0 / (_y1 - _y2);
        double db_left = (GetBValue(_color2) - GetBValue(_color0)) * 1.0 / (_y2 - _y0);
        double db_right = (GetBValue(_color1) - GetBValue(_color2)) * 1.0 / (_y1 - _y2);


        //开始进行填充
        double xs = _x0;
        double xe = _x1;
        double rs = GetRValue(_color0);
        double gs = GetGValue(_color0);
        double bs = GetBValue(_color0);
        double re = GetRValue(_color1);
        double ge = GetGValue(_color1);
        double be = GetBValue(_color1);

        if (abs(_ver2.pos.y - _ver0.pos.y) < 0.0001f)
            return;

        for (double y = _y0; y <= _y2; y += 1.0f)
        {
            double s = (double)(y - _ver0.pos.y) /
                (_ver2.pos.y - _ver0.pos.y);
            double Zt = 0.0;
            double t = s;
            double Z0 = _ver0.zdepth;
            double Z2 = _ver2.zdepth;

            if (Scene::GetInstance()->CorrectFlag) {
                if (Z0 != 0.0 && Z2 != 0.0) Zt = 1 / Z0 + s * (1 / Z2 - 1 / Z0);
                if (Zt != 0) Zt = 1 / Zt;
                if (Z2 != Z0) t = (Zt - Z0) / (Z2 - Z0);
            }

            double lz = Zt;
            double lu = _ver0.tc.u +
                t * (_ver2.tc.u - _ver0.tc.u);
            double lv = _ver0.tc.v +
                t * (_ver2.tc.v - _ver0.tc.v);

            t = s;
            double Z1 = _ver1.zdepth;
            if (Scene::GetInstance()->CorrectFlag) {
                if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
                if (Zt != 0) Zt = 1 / Zt;
                if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
            }

            double rz = Zt;
            double ru = _ver1.tc.u +
                t * (_ver2.tc.u - _ver1.tc.u);
            double rv = _ver1.tc.v +
                t * (_ver2.tc.v - _ver1.tc.v);

            DDALine(xs, xe, y,
                rs, gs, bs, re, ge, be,
                lz, lu, lv, rz, ru, rv);
            xs += dxy_left;
            xe += dxy_right;
            rs += dr_left;
            gs += dg_left;
            bs += db_left;
            re += dr_right;
            ge += dg_right;
            be += db_right;
        }
    }

    void DrawBottomTriangle(Vertex& _ver0, Vertex& _ver1, Vertex& _ver2) {
        //判断输入的三角形
        if (abs(_ver2.pos.y - _ver1.pos.y) < 0.000001f)
        {/*empty*/
        }
        else if (abs(_ver0.pos.y - _ver2.pos.y) < 0.000001f)
        {
            SwapVertex(_ver0, _ver1);
        }
        else if (abs(_ver1.pos.y - _ver0.pos.y) < 0.000001f)
        {
            SwapVertex(_ver0, _ver2);
        }
        else
        { //不是平底三角形
            return;
        }
        if (_ver1.pos.x < _ver2.pos.x)
        {
            SwapVertex(_ver1, _ver2);
        }
        else if (abs(_ver1.pos.x - _ver2.pos.x) < 0.000001f)
        { //不是三角形
            return;
        }

        double _x0 = _ver0.pos.x;
        double _x1 = _ver1.pos.x;
        double _x2 = _ver2.pos.x;
        double _y0 = _ver0.pos.y;
        double _y1 = _ver1.pos.y;
        double _y2 = _ver2.pos.y;
        COLORREF _color0 = _ver0.color;
        COLORREF _color1 = _ver1.color;
        COLORREF _color2 = _ver2.color;

        //计算左右误差
        double dxy_left = (_x2 - _x0) * 1.0 / (_y2 - _y0);
        double dxy_right = (_x1 - _x0) * 1.0 / (_y1 - _y0);
        double dr_left = (GetRValue(_color2) - GetRValue(_color0)) * 1.0 / (_y2 - _y0);
        double dr_right = (GetRValue(_color1) - GetRValue(_color0)) * 1.0 / (_y1 - _y0);
        double dg_left = (GetGValue(_color2) - GetGValue(_color0)) * 1.0 / (_y2 - _y0);
        double dg_right = (GetGValue(_color1) - GetGValue(_color0)) * 1.0 / (_y1 - _y0);
        double db_left = (GetBValue(_color2) - GetBValue(_color0)) * 1.0 / (_y2 - _y0);
        double db_right = (GetBValue(_color1) - GetBValue(_color0)) * 1.0 / (_y1 - _y0);


        //开始进行填充
        double xs = _x0;
        double xe = _x0;
        double rs = GetRValue(_color0);
        double gs = GetGValue(_color0);
        double bs = GetBValue(_color0);
        double re = GetRValue(_color0);
        double ge = GetGValue(_color0);
        double be = GetBValue(_color0);

        if (abs(_ver2.pos.y - _ver0.pos.y) < 0.0001f)
            return;

        for (double y = _y0; y <= _y2; y += 1.0f)
        {
            double s = (double)(y - _ver0.pos.y) /
                (_ver2.pos.y - _ver0.pos.y);
            double Zt = 0.0;
            double t = s;
            double Z0 = _ver0.zdepth;
            double Z2 = _ver2.zdepth;
            
            if (Scene::GetInstance()->CorrectFlag) {
                if (Z0 != 0.0 && Z2 != 0.0) Zt = 1 / Z0 + s * (1 / Z2 - 1 / Z0);
                if (Zt != 0) Zt = 1 / Zt;
                if (Z2 != Z0) t = (Zt - Z0) / (Z2 - Z0);
            }

            double lz = Zt;
            double lu = _ver0.tc.u +
                t * (_ver2.tc.u - _ver0.tc.u);
            double lv = _ver0.tc.v +
                t * (_ver2.tc.v - _ver0.tc.v);

            t = s;
            double Z1 = _ver1.zdepth;
            if (Scene::GetInstance()->CorrectFlag) {
                if (Z0 != 0.0 && Z1 != 0.0) Zt = 1 / Z0 + s * (1 / Z1 - 1 / Z0);
                if (Zt != 0) Zt = 1 / Zt;
                if (Z1 != Z0) t = (Zt - Z0) / (Z1 - Z0);
            }

            double rz = Zt;
            double ru = _ver0.tc.u +
                t * (_ver1.tc.u - _ver0.tc.u);
            double rv = _ver0.tc.v +
                t * (_ver1.tc.v - _ver0.tc.v);

            DDALine(xs, xe, y,
                rs, gs, bs, re, ge, be,
                lz, lu, lv, rz, ru, rv);
            xs += dxy_left;
            xe += dxy_right;
            rs += dr_left;
            gs += dg_left;
            bs += db_left;
            re += dr_right;
            ge += dg_right;
            be += db_right;
        }
    }

    void Draw2DTriangle(Vertex& _ver0, Vertex& _ver1, Vertex& _ver2) {
        int _x0 = _ver0.pos.x;
        int _x1 = _ver1.pos.x;
        int _x2 = _ver2.pos.x;
        int _y0 = _ver0.pos.y;
        int _y1 = _ver1.pos.y;
        int _y2 = _ver2.pos.y;

        if ((_x0 == _x1 && _x1 == _x2)
            || (_y0 == _y1 && _y1 == _y2))
        {
            return; //传入的点无法构成三角形
        }

        //将三个点按自上而下排序
        if (_ver0.pos.y > _ver1.pos.y)
        {
            SwapVertex(_ver0, _ver1);
        }
        if (_ver0.pos.y > _ver2.pos.y)
        {
            SwapVertex(_ver0, _ver2);
        }
        if (_ver1.pos.y > _ver2.pos.y)
        {
            SwapVertex(_ver1, _ver2);
        }

        //进行绘制
        if (abs(_ver0.pos.y - _ver1.pos.y) < 0.000001f) //平顶三角形
        {
            DrawTopTriangle(_ver0, _ver1, _ver2);
        }
        else if (abs(_ver1.pos.y - _ver2.pos.y) < 0.000001f) //平底三角形
        {
            DrawBottomTriangle(_ver0, _ver1, _ver2);
        }
        else
        { //分成一个平底三角形和一个平顶三角形
            int newX = _ver0.pos.x + 0.5 +
                (double)1.0 * (_ver1.pos.y - _ver0.pos.y) *
                (_ver2.pos.x - _ver0.pos.x) / (_ver2.pos.y - _ver0.pos.y);
            double newR = GetRValue(_ver0.color) +
                (double)1.0 * (_ver1.pos.y - _ver0.pos.y) *
                (GetRValue(_ver2.color) - GetRValue(_ver0.color)) /
                (_ver2.pos.y - _ver0.pos.y);
            double newG = GetGValue(_ver0.color) +
                (double)1.0 * (_ver1.pos.y - _ver0.pos.y) *
                (GetGValue(_ver2.color) - GetGValue(_ver0.color)) /
                (_ver2.pos.y - _ver0.pos.y);
            double newB = GetBValue(_ver0.color) +
                (double)1.0 * (_ver1.pos.y - _ver0.pos.y) *
                (GetBValue(_ver2.color) - GetBValue(_ver0.color)) /
                (_ver2.pos.y - _ver0.pos.y);
            //新点材质的颜色 线性插值 未做透视矫正
            COLORREF newColor = RGB((int)newR, (int)newG, (int)newB);
            Vertex newVertex = _ver1;
            newVertex.pos.x = newX;
            newVertex.color = newColor;
            //todo: 插值得到新点的深度和纹理坐标 做透视矫正

            double s = (_ver1.pos.y - _ver0.pos.y) /
                (_ver2.pos.y - _ver0.pos.y);

            double Z1 = 0.0;
            double Z2 = 0.0;
            double Zt = 0.0;
            double t = s;
            Z1 = _ver0.zdepth;
            Z2 = _ver2.zdepth;

            if (Scene::GetInstance()->CorrectFlag) {
                if (Z1 != 0.0 && Z2 != 0.0) Zt = 1 / Z1 + s * (1 / Z2 - 1 / Z1);
                if (Zt != 0) Zt = 1 / Zt;
                if (Z2 != Z1) t = (Zt - Z1) / (Z2 - Z1);
            }

            newVertex.zdepth = Zt;
            newVertex.tc.u = _ver0.tc.u +
                t * (_ver2.tc.u - _ver0.tc.u);
            newVertex.tc.v = _ver0.tc.v +
                t * (_ver2.tc.v - _ver0.tc.v);

            DrawBottomTriangle(_ver0, newVertex, _ver1);
            DrawTopTriangle(newVertex, _ver1, _ver2);
        }
    }

    void Draw(TransformMatrix& WT, Scene* mScene) {
       // Scene* mScene = Scene::GetInstance();
        Camera& mCam = mScene->mCamera;

        TransformMatrix mObTrans = mCam.mObTrans;
        TransformMatrix mNparTrans = mCam.mNparTrans;

        Vertex temp[3];
        //mNparTrans.printTrans();
        //mObTrans.printTrans();
        //WT.printTrans();
        TransformMatrix Trans1;
        Trans1 = mNparTrans * mObTrans;
        Trans1 = Trans1 * WT;
        //Trans1.printTrans();
        for (int i = 0; i < 3; i++) {
            temp[i] = pos[i];
            temp[i].pos = Trans1 * temp[i].pos;
            temp[i].normal = Trans1 * temp[i].normal;
            temp[i].color = LightModel(temp[i]);
            temp[i].zdepth = temp[i].pos.z;
        }

        sTriangle inTri;
        vector<sTriangle> outTris;

        for (int i = 0; i < 3; i++) {
            CalculateCode(temp[i]);
            inTri.ver[i] = temp[i];
        }
        Cut(inTri, outTris);
        
        TransformMatrix mViewTrans;
        mViewTrans = mViewTrans.ViewTrans(mScene->height, mScene->width,
            mScene->maxZ, mScene->minZ, mScene->x, mScene->y);

        TransformMatrix Trans2;
        Trans2 = mViewTrans* TransformMatrix::ProjectTrans();
        //Trans2.printTrans();
        for (auto iter = outTris.begin(); iter != outTris.end(); iter++) {
            for (int i = 0; i < 3; i++) {
                iter->ver[i].pos = Trans2 * iter->ver[i].pos;
                iter->ver[i].normal = Trans2 * iter->ver[i].normal;
                //cout << iter->ver[i].pos.x << endl << iter->ver[i].pos.y << endl << iter->ver[i].pos.z << endl;
            }
            Draw2DTriangle(iter->ver[0], iter->ver[1], iter->ver[2]);
        }
    }
};

class Model {
public:
    vector<Triangle> Tris;
    TransformMatrix mWorldTrans;

    void AddTriangle(const Triangle& tri) {
        Tris.push_back(tri);
    }

    void SetWorldTrans(const TransformMatrix& WT) {
        mWorldTrans = WT;
    }

    void Draw(Scene* mScene) {
        mPixels.clear();

        for (auto iter = Tris.begin(); iter != Tris.end(); iter++) {
            iter->Draw(mWorldTrans, mScene);
        }

        HDC mHDC = GetImageHDC();
        HBITMAP bmp = CreateCompatibleBitmap(mHDC, SCREEN_Width, SCREEN_Height);
        HDC memDC = CreateCompatibleDC(mHDC);
        SelectObject(memDC, bmp);

        BITMAPINFO bmpInfo;
        bmpInfo.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
        bmpInfo.bmiHeader.biWidth = SCREEN_Width;
        bmpInfo.bmiHeader.biHeight = SCREEN_Height;
        bmpInfo.bmiHeader.biPlanes = 1;
        bmpInfo.bmiHeader.biBitCount = 24;
        bmpInfo.bmiHeader.biCompression = BI_RGB;
        bmpInfo.bmiHeader.biSizeImage = 0;
        bmpInfo.bmiHeader.biXPelsPerMeter = 3000;
        bmpInfo.bmiHeader.biYPelsPerMeter = 3000;
        bmpInfo.bmiHeader.biClrUsed = 0;
        bmpInfo.bmiHeader.biClrImportant = 0;

        BYTE red, green, blue;

        memset(DataComp, 0, mLen);
        

        for (int i = 0; i < mLen; i += 3)
        {
            DataComp[i + 2] = 0;
            DataComp[i] = 0;
            DataComp[i + 1] = 0;
        }
        for (auto iter = mPixels.begin();
            iter != mPixels.end(); ++iter)
        {
            if (iter->x >= SCREEN_Width || iter->y >= SCREEN_Height || iter->x < 0 || iter->y < 0)
                continue;
            //iter是第n个像素点 行优先 todo: 这里为什么要-y
            int n = (SCREEN_Height - iter->y) * SCREEN_Width + iter->x;
            int beginPos = n * 3; //mData中颜色开始的位置
            COLORREF color = iter->color; //当前点的颜色
            BYTE* pByte = (BYTE*)&color;
            //第一个字节为0 其它三个字节分别表示蓝色、绿色和红色
            red = *(pByte + 2);
            blue = *(pByte + 0);
            green = *(pByte + 1);
            //bmp图像中是R、G、B排列
            DataComp[beginPos] = red;
            DataComp[beginPos + 1] = green;
            DataComp[beginPos + 2] = blue;
        }
        SetDIBits(mHDC, bmp, 0, SCREEN_Height,
            DataDraw, &bmpInfo, DIB_RGB_COLORS);
        BitBlt(mHDC, 0, 0, SCREEN_Width, SCREEN_Height, memDC, 0, 0, SRCCOPY); /*将位图内容传送到目标DC上*/
        DeleteDC(memDC); /*删除内存DC*/
        DeleteObject(bmp);

        DataDraw = DataComp;

        COLORREF color = RGB(255, 255, 255);
        SetTextColor(mHDC, color);
        // int fps to LPCWSTR

        wsprintf(buffer3, L"FPS: %d", mScene->fps);
        SetBkMode(mHDC, TRANSPARENT);
        TextOut(mHDC, 5, 5, buffer0, wcslen(buffer0));
        TextOut(mHDC, 5, 25, buffer1, wcslen(buffer1));
        TextOut(mHDC, 5, 45, buffer2, wcslen(buffer2));
        TextOut(mHDC, 5, 65, buffer3, wcslen(buffer3));
    }
};

void DrawFps(int fps) {
    string str = to_string(fps);
    LPCTSTR str0 = _T("zoom\nFPS:");
    //outtextxy(5, 65, str0);
    for (int i = 0; i < str.size(); i++) {
        //outtextxy(50 + i * 8, 65, str[0]);
    }
}

int main()
{
    DataDraw = new BYTE[mLen];
    DataComp = new BYTE[mLen];
    memset(DataDraw, 0, mLen);
    memset(DataComp, 0, mLen);
    Scene* pScene = Scene::GetInstance();

    wsprintf(buffer0, L"滚轮缩放 鼠标点击拖动 操控摄像机旋转", pScene->fps);
    wsprintf(buffer1, L"W A S D 操控摄像机平面内上下左右移动", pScene->fps);
    wsprintf(buffer2, L"H 开启/关闭纹理透视修正  L 开启/关闭点光源  P 开启/关闭平行光", pScene->fps);

    Vector VRP(0, 0, 0, 1);
    Vector VPN(1, 0, 0, 1);
    Vector VUP(0 , 1, 0, 1);
    Vector PRP(5, 5, -5, 0);

    Camera cam(VPN, VRP, VUP, PRP,
        -2, 2, -2, 2, 2, 2);

    pScene->mCamera = cam;

    pScene->x = 0;
    pScene->y = 0;
    pScene->height = SCREEN_Height;
    pScene->width = SCREEN_Width;
    pScene->maxZ = 1.0;
    pScene->minZ = 0.0;

    AmbLight AL(255, 255, 255);
    DirLight DL(Vector(0, 0, 1, 0), 255, 255, 255);

    pScene->mAmbLight = AL;
    pScene->mDirLight = DL;

	Vector norm(0, 0, -1, 0);
	Material mat(0.2, 0.4, 0.6, 2);
	TexCoord tc0(1, 0), tc1(0,1), tc2(0,0), tc3(0.5, 1.7320508 / 2);

	Triangle triangle1( //三角形在模型坐标系中
		Vertex(Vector(1,0,0,1), //位置
			norm, //法向
			mat, //材质
			tc1), //纹理坐标
		Vertex(Vector(0,1,0,1),
			norm,
			mat,
			tc0), //纹理坐标
		Vertex(Vector(0,0,0,1),
			norm,
			mat,
			tc2)); //纹理坐标

    Triangle triangle2(
        Vertex(Vector(0,1,0,1),
            Vector(1,1,1,0).Normalize(),
            mat,
            tc2),
        Vertex(Vector(1, 0, 0, 1),
            Vector(1, 1, 1, 0).Normalize(),
            mat,
            tc3),
        Vertex(Vector(0, 0, 1, 1),
            Vector(1, 1, 1, 0).Normalize(),
            mat,
            tc0));

    Triangle triangle3(
        Vertex(Vector(0, 0, 1, 1),
            Vector(-1,0,0,0),
            mat,
            tc1),
        Vertex(Vector(0, 0, 0, 1),
            Vector(-1, 0, 0, 0),
            mat,
            tc2),
        Vertex(Vector(0, 1, 0, 1),
            Vector(-1, 0, 0, 0),
            mat,
            tc0));

    Triangle triangle4(
        Vertex(Vector(1, 0, 0, 1),
            Vector(0, -1, 0, 0),
            mat,
            tc1),
        Vertex(Vector(0, 0, 0, 1),
            Vector(0, -1, 0, 0),
            mat,
            tc2),
        Vertex(Vector(0, 0, 1, 1),
            Vector(0, -1, 0, 0),
            mat,
            tc0));

    Model mModel;
    mModel.AddTriangle(triangle1);
    mModel.AddTriangle(triangle2);
    mModel.AddTriangle(triangle3);
    mModel.AddTriangle(triangle4);

    Vector modelPos(0, 0, 0, 1);
    Vector scale(1, 1, 1, 1);

    TransformMatrix WT;
    WT = WT.WorldTrans(modelPos, scale, 0.0f, 0.0f, 0.0f);
    mModel.mWorldTrans = WT;

    pScene->ReadTexture();
    Scene::mInstance = pScene;

    initgraph(1000, 800);
    settextcolor(WHITE);
    settextstyle(20, 0, _T("Consolas"));
    setlocale(LC_ALL, "chs");

    int frame = 0;

    ExMessage msg;

    while (1) {
        time_t lastTime = time(NULL);

        mModel.Draw(pScene);

        msg = getmessage(EM_KEY | EM_MOUSE);
        if (msg.message == WM_KEYDOWN)
        {
            switch (msg.vkcode)
            {
            case 37://left
            {
                pScene->mCamera.mObTrans = WT.RotateAxisTrans(VPN.Normalize(), 5) * pScene->mCamera.mObTrans;
                ///mModel.mWorldTrans = WT;
               // cleardevice();
                break;
            }
            case 38://up
            {
                Vector UAxis = CrossProduct(VUP, VPN).Normalize();
                Vector VAxis = CrossProduct(VPN.Normalize(), UAxis);
                pScene->mCamera.mObTrans = WT.RotateAxisTrans(VAxis, 5) * pScene->mCamera.mObTrans;
                //mModel.mWorldTrans = WT;
               // cleardevice();
                break;
            }
            case 39://down
            {
                pScene->mCamera.mObTrans = WT.RotateZTrans(-10) * pScene->mCamera.mObTrans;
                //mModel.mWorldTrans = WT;
               // cleardevice();
                break;
            }
            case 40://right
            {
                pScene->mCamera.mObTrans = WT.RotateAxisTrans(CrossProduct(VUP, VPN).Normalize(), 5) * pScene->mCamera.mObTrans;
                //mModel.mWorldTrans = WT;
                //cleardevice();
                break;
            }
            case 65://D
            {
                pScene->mCamera.mObTrans = WT.MoveTrans(Vector(0.05, 0, 0)) * pScene->mCamera.mObTrans;
              //  cleardevice();
                break;
            }
            case 68://A
            {
                pScene->mCamera.mObTrans = WT.MoveTrans(Vector(-0.05, 0, 0)) * pScene->mCamera.mObTrans;
               // cleardevice();
                break;
            }
            case 87://W
            {
                pScene->mCamera.mObTrans = WT.MoveTrans(Vector(0, -0.05, 0)) * pScene->mCamera.mObTrans;
              //  cleardevice();
                break;
            }
            case 83://S
            {
                pScene->mCamera.mObTrans = WT.MoveTrans(Vector(0, 0.05, 0)) * pScene->mCamera.mObTrans;
              //  cleardevice();
                break;
            }
            case 0x48: //H
            {
                pScene->CorrectFlag = 1 - pScene->CorrectFlag;
               // cleardevice();
                break;
            }
            case 0x4C: //L
            {
                AL = AmbLight(255-AL.r, 255-AL.g, 255-AL.b);
                pScene->mAmbLight = AL;
                //cleardevice();
                break;
            }
            case 0x50: //P
            {
                DL = DirLight(Vector(0, 0, 1, 0), 255 - DL.r, 255 - DL.g, 255 - DL.b);
                pScene->mDirLight = DL;
                //cleardevice();
                break;
            }
            default:
                break;
            }
        }
        else if (EM_MOUSE) {
            switch (msg.message)
            {
            case WM_MOUSEWHEEL:
            {
                if (msg.wheel < 0) {
                    pScene->mCamera.mObTrans = WT.ZoomTrans(Vector(0.85, 0.85, 0.85, 1)) * pScene->mCamera.mObTrans;
                    //cleardevice();
                    break;
                }
                else {
                    pScene->mCamera.mObTrans = WT.ZoomTrans(Vector(1.2, 1.2, 1.2, 1)) * pScene->mCamera.mObTrans;
                    //cleardevice();
                    break;
                }
                
            }
            case WM_LBUTTONDOWN:
            {
                int x0 = msg.x;
                int y0 = msg.y;

                while (msg.message != WM_LBUTTONUP) {

                    mModel.Draw(pScene);
                    msg = getmessage(EM_MOUSE);
                    if (msg.message == WM_MOUSEMOVE) {
                        int x1 = msg.x;
                        int y1 = msg.y;
                        Vector move(y0 - y1, x0 - x1, 0, 0);
                        pScene->mCamera.mObTrans = WT.RotateAxisTrans(move.Reverse().Normalize(), move.length() * 0.02) * pScene->mCamera.mObTrans;
                        //cleardevice();
                    }
                }
                break;
            }
            default:
                break;
            }
        }

        ++frame;
        time_t deltaTime = time(NULL) - lastTime;
        if (deltaTime > 0) {
            lastTime = time(NULL);
            pScene->fps = frame;
            frame = 0;
        }
       
    }

    closegraph();
    return 0;
}

