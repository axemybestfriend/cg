#include <GL/glut.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

#define M_PI 3.14159265f

float Ka = 0.2f;
float Kd = 0.7f;
float Ks = 0.2f;
float n = 8.0f;
float Ia = 0.3f;
float Il = 5.0f;
float k = 1.0f;


struct Point {
    float x, y, z;
    Point(float x = 0, float y = 0, float z = 0) : x(x), y(y), z(z) {}

    Point operator+(const Point& o) const { return Point(x + o.x, y + o.y, z + o.z); }
    Point operator-(const Point& o) const { return Point(x - o.x, y - o.y, z - o.z); }
    Point operator*(float s) const { return Point(x * s, y * s, z * s); }
};

Point lightPos(0, 5, 5);
Point viewPos(0, 2, 6);

struct Triangle {
    int v1, v2, v3; 
    Triangle(int v1, int v2, int v3) : v1(v1), v2(v2), v3(v3) {} 
};

std::vector<Point> vertices, vertexNormals;
std::vector<float> vertexIntensities;
std::vector<Triangle> triangles;

float distance(const Point& a, const Point& b) {
    float dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

Point normalize(const Point& v) {
    float len = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    if (len <= 0.0001f) return Point(0, 1, 0);
    return (len > 0) ? Point(v.x / len, v.y / len, v.z / len) : v;
}

float scalar(const Point& a, const Point& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

Point cross(const Point& a, const Point& b) {
    return Point(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

int profileSegments = 18;
int rotationSegments = 72;
float radius = 1.0f;
std::vector<Point> profile;

void createFigure() {
    profile.clear();
    for (int i = 0; i <= profileSegments; i++) {
        float t = M_PI * i / profileSegments;
        profile.push_back(Point(radius * sin(t), radius * cos(t), 0));
    }
}

std::vector<Point> rotateProfile(const std::vector<Point>& prof, float a) {
    std::vector<Point> r;
    float ca = cos(a), sa = sin(a);
    for (const auto& p : prof) {
        float x = p.x * ca - p.z * sa;
        float z = p.x * sa + p.z * ca;
        r.push_back(Point(x, p.y, z));
    }
    return r;
}

Point triNormal(const Point& a, const Point& b, const Point& c) {
    return normalize(cross(b - a, c - a));
}

float calcIntensity(const Point& normal, const Point& point) {
    Point L = normalize(lightPos - point);
    Point S = normalize(viewPos - point);

    float NL = std::max(0.0f, scalar(normal, L));

    Point R = normalize(normal * (2.0f * NL) - L);
    float RS = std::max(0.0f, scalar(R, S));

    float d = distance(point, lightPos);
    return Ia * Ka + (Il / (k + d)) * (Kd * NL + Ks * pow(RS, n));
}

void buildMesh() {
    vertices.clear();
    triangles.clear();

    std::vector<int> prev;

    for (int r = 0; r <= rotationSegments; r++) {
        float ang = 2 * M_PI * r / rotationSegments;
        auto ring = rotateProfile(profile, ang);
        std::vector<int> curr;

        for (const auto& p : ring) {
            vertices.push_back(p);
            curr.push_back(vertices.size() - 1);
        }

        if (!prev.empty()) {
            for (int i = 0; i < profileSegments; i++) {
                int a = prev[i], b = prev[i + 1], c = curr[i + 1], d = curr[i];
                triangles.emplace_back(a, b, c);
                triangles.emplace_back(a, c, d);
            }
        }
        prev = curr;
    }

    vertexNormals.assign(vertices.size(), Point(0, 0, 0));

    for (auto& tri : triangles) {
        Point a = vertices[tri.v1];
        Point b = vertices[tri.v2];
        Point c = vertices[tri.v3];

        Point triNormal = normalize(cross(c - a, b - a));

        float length = sqrt(triNormal.x * triNormal.x +
            triNormal.y * triNormal.y +
            triNormal.z * triNormal.z);

        if (length > 0.0001f) {
            triNormal = triNormal * (1.0f / length);

            vertexNormals[tri.v1] = vertexNormals[tri.v1] + triNormal;
            vertexNormals[tri.v2] = vertexNormals[tri.v2] + triNormal;
            vertexNormals[tri.v3] = vertexNormals[tri.v3] + triNormal;
        }
    }

    for (size_t i = 0; i < vertices.size(); i++) {
        vertexNormals[i] = normalize(vertexNormals[i]);
    }

    vertexIntensities.resize(vertices.size());
    for (size_t i = 0; i < vertices.size(); i++)
        vertexIntensities[i] = calcIntensity(vertexNormals[i], vertices[i]);
}

float yaw = -90.0f;
float pitch = 0.0f;
Point front(0, 0, -1);
Point worldUp(0, 1, 0);
float speed = 0.2f;

Point right, up;

void updateCameraVectors() {
    front.x = cosf(yaw * M_PI / 180.0f) * cosf(pitch * M_PI / 180.0f);
    front.y = sinf(pitch * M_PI / 180.0f);
    front.z = sinf(yaw * M_PI / 180.0f) * cosf(pitch * M_PI / 180.0f);
    front = normalize(front);

    Point horizontalFront = front;
    horizontalFront.y = 0;
    right = normalize(cross(horizontalFront, worldUp));

    up = normalize(cross(right, front));
}


void keyboard(unsigned char key, int x, int y) {
    if (key == 'w') viewPos = viewPos + front * speed;
    if (key == 's') viewPos = viewPos - front * speed;
    if (key == 'a') viewPos = viewPos - right * speed;
    if (key == 'd') viewPos = viewPos + right * speed;

    glutPostRedisplay();
}


void drawLineInterpolated(const Point& start, const Point& end,
    float iStart, float iEnd) {
    int steps = 100;
    for (int i = 0; i <= steps; i++) {
        float t = (float)i / steps;
        Point p = start + (end - start) * t;
        float intensity = iStart + (iEnd - iStart) * t;
        intensity = std::max(0.0f, std::min(1.0f, intensity));

        glColor3f(intensity, intensity, intensity);
        glVertex3f(p.x, p.y, p.z);
    }
}

int lastMouseX = 400, lastMouseY = 300;
bool firstMouse = true;

void mouseMotion(int x, int y) {
    int centerX = 400, centerY = 300;
    static bool firstMouse = true;

    if (firstMouse) {
        lastMouseX = x;
        lastMouseY = y;
        firstMouse = false;
    }

    float xoffset = x - lastMouseX;
    float yoffset = lastMouseY - y;
    lastMouseX = x;
    lastMouseY = y;

    float sensitivity = 0.1f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    yaw += xoffset;
    pitch += yoffset;

    if (pitch > 89.0f) pitch = 89.0f;
    if (pitch < -89.0f) pitch = -89.0f;

    updateCameraVectors();

    glutWarpPointer(centerX, centerY);
    lastMouseX = centerX;
    lastMouseY = centerY;

    glutPostRedisplay();
}

void drawUpperTriangle(const Point& top1, const Point& top2, const Point& bottom,
    float i_top1, float i_top2, float i_bottom) {

    // Сортируем верхние точки по X
    Point left_top, right_top;
    float i_left_top, i_right_top;

    if (top1.x <= top2.x) {
        left_top = top1; right_top = top2;
        i_left_top = i_top1; i_right_top = i_top2;
    }
    else {
        left_top = top2; right_top = top1;
        i_left_top = i_top2; i_right_top = i_top1;
    }

    // Убедимся, что bottom действительно ниже
    if (bottom.y <= left_top.y) return;

    int scanlines = (int)((bottom.y - left_top.y) / 0.01f);
    if (scanlines < 1) scanlines = 1;

    for (int i = 0; i <= scanlines; i++) {
        float y = left_top.y + (bottom.y - left_top.y) * (float)i / scanlines;

        // Интерполируем левую точку (от left_top до bottom)
        float t_left = (y - left_top.y) / (bottom.y - left_top.y);
        Point left_point = left_top + (bottom - left_top) * t_left;
        float i_left = i_left_top + (i_bottom - i_left_top) * t_left;

        // Интерполируем правую точку (от right_top до bottom)
        float t_right = (y - right_top.y) / (bottom.y - right_top.y);
        Point right_point = right_top + (bottom - right_top) * t_right;
        float i_right = i_right_top + (i_bottom - i_right_top) * t_right;

        // Рисуем горизонтальную линию между левой и правой точкой
        int points = (int)(distance(left_point, right_point) / 0.01f);
        if (points < 1) points = 1;

        for (int j = 0; j <= points; j++) {
            float t = (float)j / points;
            Point p = left_point + (right_point - left_point) * t;
            float intensity = i_left + (i_right - i_left) * t;
            intensity = std::max(0.0f, std::min(1.0f, intensity));

            glColor3f(intensity, intensity, intensity);
            glVertex3f(p.x, p.y, p.z);
        }
    }
}

void drawLowerTriangle(const Point& top, const Point& bottom1, const Point& bottom2,
    float i_top, float i_bottom1, float i_bottom2) {

    Point left_bottom, right_bottom;
    float i_left_bottom, i_right_bottom;

    if (bottom1.x <= bottom2.x) {
        left_bottom = bottom1; right_bottom = bottom2;
        i_left_bottom = i_bottom1; i_right_bottom = i_bottom2;
    }
    else {
        left_bottom = bottom2; right_bottom = bottom1;
        i_left_bottom = i_bottom2; i_right_bottom = i_bottom1;
    }

    // Убедимся, что top действительно выше
    if (top.y >= left_bottom.y) return;

    int scanlines = (int)((left_bottom.y - top.y) / 0.01f);
    if (scanlines < 1) scanlines = 1;

    for (int i = 0; i <= scanlines; i++) {
        float y = top.y + (left_bottom.y - top.y) * (float)i / scanlines;

        // Интерполируем левую точку (от top до left_bottom)
        float t_left = (y - top.y) / (left_bottom.y - top.y);
        Point left_point = top + (left_bottom - top) * t_left;
        float i_left = i_top + (i_left_bottom - i_top) * t_left;

        // Интерполируем правую точку (от top до right_bottom)
        float t_right = (y - top.y) / (right_bottom.y - top.y);
        Point right_point = top + (right_bottom - top) * t_right;
        float i_right = i_top + (i_right_bottom - i_top) * t_right;

        // Рисуем горизонтальную линию между левой и правой точкой
        int points = (int)(distance(left_point, right_point) / 0.01f);
        if (points < 1) points = 1;

        for (int j = 0; j <= points; j++) {
            float t = (float)j / points;
            Point p = left_point + (right_point - left_point) * t;
            float intensity = i_left + (i_right - i_left) * t;
            intensity = std::max(0.0f, std::min(1.0f, intensity));

            glColor3f(intensity, intensity, intensity);
            glVertex3f(p.x, p.y, p.z);
        }
    }
}

void drawTriangleGouraud(const Point& v1, const Point& v2, const Point& v3,
    float i1, float i2, float i3) {

    // Сначала рисуем границы для отладки
    drawLineInterpolated(v1, v2, i1, i2);
    drawLineInterpolated(v2, v3, i2, i3);
    drawLineInterpolated(v3, v1, i3, i1);

    // Определяем тип треугольника
    float epsilon = 0.001f;

    if (fabs(v1.y - v2.y) < epsilon) {
        // Верхний треугольник: v1 и v2 сверху
        if (v3.y > v1.y) { // v3 должен быть ниже
            drawUpperTriangle(v1, v2, v3, i1, i2, i3);
        }
    }
    else if (fabs(v2.y - v3.y) < epsilon) {
        // Нижний треугольник: v2 и v3 снизу
        if (v1.y < v2.y) { // v1 должен быть выше
            drawLowerTriangle(v1, v2, v3, i1, i2, i3);
        }
    }
    else if (fabs(v1.y - v3.y) < epsilon) {
        // Две точки с одинаковым Y
        if (v1.y > v2.y) {
            // v1 и v3 сверху, v2 снизу
            drawUpperTriangle(v1, v3, v2, i1, i3, i2);
        }
        else {
            // v1 и v3 снизу, v2 сверху
            drawLowerTriangle(v2, v1, v3, i2, i1, i3);
        }
    }
    else {
        // Общий случай - сортируем по Y и разбиваем
        std::vector<std::pair<Point, float>> points = { {v1, i1}, {v2, i2}, {v3, i3} };
        std::sort(points.begin(), points.end(),
            [](const auto& a, const auto& b) { return a.first.y < b.first.y; });

        Point top = points[0].first;
        Point mid = points[1].first;
        Point bot = points[2].first;
        float i_top = points[0].second;
        float i_mid = points[1].second;
        float i_bot = points[2].second;

        // Точка разделения на стороне top-bot
        float t = (mid.y - top.y) / (bot.y - top.y);
        Point split = top + (bot - top) * t;
        float i_split = i_top + (i_bot - i_top) * t;

        // Верхний треугольник
        drawUpperTriangle(top, mid, split, i_top, i_mid, i_split);
        // Нижний треугольник
        drawLowerTriangle(mid, split, bot, i_mid, i_split, i_bot);
    }
}
void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    Point center = viewPos + front;
    gluLookAt(viewPos.x, viewPos.y, viewPos.z, center.x, center.y, center.z, up.x, up.y, up.z);
   
    updateCameraVectors();
    glPointSize(2.0f);

    glBegin(GL_TRIANGLES);
    for (auto& t : triangles) {
        glColor3f(vertexIntensities[t.v1], vertexIntensities[t.v1], vertexIntensities[t.v1]);
        glVertex3f(vertices[t.v1].x, vertices[t.v1].y, vertices[t.v1].z);

        glColor3f(vertexIntensities[t.v2], vertexIntensities[t.v2], vertexIntensities[t.v2]);
        glVertex3f(vertices[t.v2].x, vertices[t.v2].y, vertices[t.v2].z);

        glColor3f(vertexIntensities[t.v3], vertexIntensities[t.v3], vertexIntensities[t.v3]);
        glVertex3f(vertices[t.v3].x, vertices[t.v3].y, vertices[t.v3].z);
    }
    glEnd();

    glutSwapBuffers();
}

void reshape(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, float(w) / h, 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1000, 1000);
    glutCreateWindow("lab3");

    glEnable(GL_DEPTH_TEST);
    glClearColor(0, 0, 0, 1);

    createFigure();
    buildMesh();


    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutPassiveMotionFunc(mouseMotion);

    glutSetCursor(GLUT_CURSOR_NONE);


    glutMainLoop();
    return 0;
}
