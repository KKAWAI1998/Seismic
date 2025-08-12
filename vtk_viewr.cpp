// Minimal VTK PolyData viewer (fast): SDL2 + OpenGL (compat) + VBO + Scalars streaming
// 操作: 左ドラッグ=回転 / ホイール=ズーム / Space=再生/停止 / ←→=前後
//       C=カラーマップ / A=オートスケール / B=背景 / W=ワイヤ / R=リセット / O=開く / Esc=終了
#define SDL_MAIN_HANDLED
#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>

#include <commdlg.h>
#endif

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h> // ← PFNGL* 型はこれで入る
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// ---------- small helpers ----------
static inline std::string up(const std::string &s) {
    std::string t = s;
    for (char &c : t)
        c = (char)std::toupper((unsigned char)c);
    return t;
}
static inline bool fileExists(const std::string &p) {
    std::ifstream f(p);
    return (bool)f;
}

// ---------- mesh ----------
struct Mesh {
    std::vector<float> pos;    // N*3
    std::vector<unsigned> tri; // 3*M
    std::vector<float> scal;   // N
    int nvert = 0, ntri = 0;
    float smin = 0.f, smax = 1.f;
    bool ok() const { return nvert > 0 && !tri.empty(); }
};

// ---------- VTK parse (1st frame: geometry+scalars) ----------
static bool parseVTK_full(const std::string &path, Mesh &out) {
    std::ifstream ifs(path);
    if (!ifs) {
        std::cerr << "Cannot open: " << path << "\n";
        return false;
    }
    std::string tok;
    int Npts = 0;
    std::vector<float> P;
    std::vector<unsigned> T;
    std::vector<float> S;
    while (ifs >> tok) {
        std::string TUP = up(tok);
        if (TUP == "POINTS") {
            std::string dtype;
            ifs >> Npts >> dtype;
            if (Npts <= 0)
                return false;
            P.resize(Npts * 3);
            for (int i = 0; i < Npts * 3; ++i)
                ifs >> P[i];
        } else if (TUP == "POLYGONS") {
            int M = 0, tot = 0;
            ifs >> M >> tot;
            for (int k = 0; k < M; ++k) {
                int n = 0;
                ifs >> n;
                std::vector<unsigned> v(n);
                for (int i = 0; i < n; ++i)
                    ifs >> v[i];
                for (int i = 1; i + 1 < n; ++i) {
                    T.push_back(v[0]);
                    T.push_back(v[i]);
                    T.push_back(v[i + 1]);
                }
            }
        } else if (TUP == "POINT_DATA") {
            int N = 0;
            ifs >> N;
            std::streampos pos = ifs.tellg();
            std::string k;
            ifs >> k;
            if (up(k) != "SCALARS") {
                ifs.seekg(pos);
                continue;
            }
            std::string name, dtype;
            int ncomp = 1;
            ifs >> name >> dtype >> ncomp;
            std::streampos pos2 = ifs.tellg();
            std::string maybe;
            ifs >> maybe;
            if (up(maybe) != "LOOKUP_TABLE")
                ifs.seekg(pos2);
            else {
                std::string dummy;
                ifs >> dummy;
            }
            S.resize(N);
            for (int i = 0; i < N; ++i)
                ifs >> S[i];
        } else if (TUP == "FIELD") {
            std::string f;
            int m = 0;
            ifs >> f >> m;
            for (int a = 0; a < m; ++a) {
                std::string nm, ty;
                int nc = 0, nt = 0;
                ifs >> nm >> nc >> nt >> ty;
                double tmp;
                for (long long i = 0; i < (long long)nc * nt; ++i)
                    ifs >> tmp;
            }
        }
    }
    if (P.empty() || T.empty()) {
        std::cerr << "VTK missing geometry\n";
        return false;
    }
    if (S.size() != (size_t)(P.size() / 3))
        S.assign(P.size() / 3, 0.f);

    float bmin[3] = {+1e30f, +1e30f, +1e30f}, bmax[3] = {-1e30f, -1e30f, -1e30f};
    for (size_t i = 0; i < P.size(); i += 3) {
        bmin[0] = std::min(bmin[0], P[i + 0]);
        bmax[0] = std::max(bmax[0], P[i + 0]);
        bmin[1] = std::min(bmin[1], P[i + 1]);
        bmax[1] = std::max(bmax[1], P[i + 1]);
        bmin[2] = std::min(bmin[2], P[i + 2]);
        bmax[2] = std::max(bmax[2], P[i + 2]);
    }
    float c[3] = {(bmin[0] + bmax[0]) * 0.5f, (bmin[1] + bmax[1]) * 0.5f,
                  (bmin[2] + bmax[2]) * 0.5f};
    float r = 0.f;
    for (size_t i = 0; i < P.size(); i += 3) {
        float dx = P[i] - c[0], dy = P[i + 1] - c[1], dz = P[i + 2] - c[2];
        r = std::max(r, std::sqrt(dx * dx + dy * dy + dz * dz));
    }
    if (r <= 0)
        r = 1.f;
    out.pos.resize(P.size());
    for (size_t i = 0; i < P.size(); i += 3) {
        out.pos[i] = (P[i] - c[0]) / r;
        out.pos[i + 1] = (P[i + 1] - c[1]) / r;
        out.pos[i + 2] = (P[i + 2] - c[2]) / r;
    }
    out.tri.swap(T);
    out.scal.swap(S);
    out.nvert = (int)out.pos.size() / 3;
    out.ntri = (int)out.tri.size() / 3;
    auto [mn, mx] = std::minmax_element(out.scal.begin(), out.scal.end());
    out.smin = out.scal.empty() ? 0.f : *mn;
    out.smax = out.scal.empty() ? 1.f : *mx;
    std::cout << "Loaded " << path << " verts=" << out.nvert << " tris=" << out.ntri << " u=["
              << out.smin << "," << out.smax << "]\n";
    return true;
}

// ---------- VTK parse (2nd~: scalars only) ----------
static bool parseVTK_scalarsOnly(const std::string &path, int expectedN, std::vector<float> &S,
                                 float &smin, float &smax) {
    std::ifstream ifs(path);
    if (!ifs)
        return false;
    std::string w;
    int N = -1;
    while (ifs >> w) {
        if (up(w) == "POINT_DATA") {
            ifs >> N;
            break;
        }
    }
    if (N <= 0 || (expectedN > 0 && N != expectedN)) {
        std::cerr << "POINT_DATA mismatch " << path << "\n";
        return false;
    }
    while (ifs >> w) {
        if (up(w) == "SCALARS")
            break;
    }
    if (!ifs)
        return false;
    std::string name, dtype;
    int ncomp = 1;
    ifs >> name >> dtype >> ncomp;
    std::streampos pos2 = ifs.tellg();
    std::string maybe;
    ifs >> maybe;
    if (up(maybe) != "LOOKUP_TABLE")
        ifs.seekg(pos2);
    else {
        std::string dummy;
        ifs >> dummy;
    }
    S.resize(N);
    smin = +1e30f;
    smax = -1e30f;
    for (int i = 0; i < N; ++i) {
        ifs >> S[i];
        smin = std::min(smin, S[i]);
        smax = std::max(smax, S[i]);
    }
    return (bool)ifs;
}

// ---------- colormap ----------
static void mapGray(float t, float &r, float &g, float &b) {
    float c = std::clamp(t, 0.f, 1.f);
    r = g = b = c;
}
static void mapRdBu(float t, float &r, float &g, float &b) {
    static const float pal[11][3] = {
        {0.0196f, 0.1882f, 0.3804f}, {0.1294f, 0.4f, 0.6745f},    {0.2627f, 0.5765f, 0.7647f},
        {0.5725f, 0.7725f, 0.8706f}, {0.8196f, 0.8980f, 0.9412f}, {0.9686f, 0.9686f, 0.9686f},
        {0.9922f, 0.8588f, 0.7804f}, {0.9882f, 0.6824f, 0.5529f}, {0.9569f, 0.4275f, 0.2627f},
        {0.8392f, 0.1529f, 0.1569f}, {0.6471f, 0.0f, 0.1490f}};
    t = std::clamp(t, 0.f, 1.f) * 10.f;
    int i = (int)t;
    float f = t - i;
    if (i >= 10) {
        i = 10;
        f = 0;
    }
    r = pal[i][0] * (1 - f) + pal[std::min(i + 1, 10)][0] * f;
    g = pal[i][1] * (1 - f) + pal[std::min(i + 1, 10)][1] * f;
    b = pal[i][2] * (1 - f) + pal[std::min(i + 1, 10)][2] * f;
}

// ---------- VBO loader ----------
static PFNGLGENBUFFERSPROC pglGenBuffers = nullptr;
static PFNGLBINDBUFFERPROC pglBindBuffer = nullptr;
static PFNGLBUFFERDATAPROC pglBufferData = nullptr;
static PFNGLBUFFERSUBDATAPROC pglBufferSubData = nullptr;
static PFNGLDELETEBUFFERSPROC pglDeleteBuffers = nullptr;
static void loadVBOFuncs() {
    pglGenBuffers = (PFNGLGENBUFFERSPROC)SDL_GL_GetProcAddress("glGenBuffers");
    pglBindBuffer = (PFNGLBINDBUFFERPROC)SDL_GL_GetProcAddress("glBindBuffer");
    pglBufferData = (PFNGLBUFFERDATAPROC)SDL_GL_GetProcAddress("glBufferData");
    pglBufferSubData = (PFNGLBUFFERSUBDATAPROC)SDL_GL_GetProcAddress("glBufferSubData");
    pglDeleteBuffers = (PFNGLDELETEBUFFERSPROC)SDL_GL_GetProcAddress("glDeleteBuffers");
    if (!pglGenBuffers || !pglBindBuffer || !pglBufferData || !pglBufferSubData ||
        !pglDeleteBuffers)
        std::cerr << "[WARN] VBO not available. Fallback to client arrays.\n";
}

// ---------- Viewer ----------
struct Viewer {
    std::vector<std::string> files;
    int cur = 0;
    Mesh mesh;
    std::vector<unsigned char> colors;
    GLuint vboPos = 0, vboCol = 0, ibo = 0;
    bool hasVBO() const { return pglGenBuffers && vboPos != 0; }

    bool playing = false, wire = false, bgWhite = true, autoPerFrame = true;
    enum { GRAY, RDBU } cmap = RDBU;
    int winW = 1200, winH = 900;
    float yaw = 0.f, pitch = 0.f, dist = 2.5f, fps = 24.f;
    float globalMin = 0.f, globalMax = 1.f;

    static bool parseIndexRun(const std::string &p, int &start, int &len, std::string &pref,
                              std::string &suf) {
        int dot = (int)p.find_last_of('.');
        int end = (dot == -1 ? (int)p.size() - 1 : dot - 1);
        int i = end;
        while (i >= 0 && std::isdigit((unsigned char)p[i]))
            --i;
        int j = i + 1;
        len = end - j + 1;
        if (len <= 0)
            return false;
        start = j;
        pref = p.substr(0, start);
        suf = p.substr(j + len);
        return true;
    }
    void buildAutoSeries(const std::string &sample) {
        files.clear();
        int st = 0, len = 0;
        std::string pref, suf;
        if (!parseIndexRun(sample, st, len, pref, suf)) {
            files.push_back(sample);
            return;
        }
        int curIdx = std::atoi(sample.substr(st, len).c_str());
        int first = 0;
        for (int k = curIdx; k >= 0; --k) {
            char buf[1024];
            std::snprintf(buf, sizeof(buf), ("%s%0" + std::to_string(len) + "d%s").c_str(),
                          pref.c_str(), k, suf.c_str());
            if (!fileExists(buf)) {
                first = k + 1;
                break;
            }
            if (k == 0)
                first = 0;
        }
        int last = 0;
        for (int k = curIdx;; ++k) {
            char buf[1024];
            std::snprintf(buf, sizeof(buf), ("%s%0" + std::to_string(len) + "d%s").c_str(),
                          pref.c_str(), k, suf.c_str());
            if (!fileExists(buf)) {
                last = k - 1;
                break;
            }
        }
        for (int k = first; k <= last; ++k) {
            char buf[1024];
            std::snprintf(buf, sizeof(buf), ("%s%0" + std::to_string(len) + "d%s").c_str(),
                          pref.c_str(), k, suf.c_str());
            files.emplace_back(buf);
        }
        if (files.empty())
            files.push_back(sample);
    }

    void setupGLBuffersOnce() {
        if (!mesh.ok() || vboPos)
            return;
        if (pglGenBuffers) {
            pglGenBuffers(1, &vboPos);
            pglBindBuffer(GL_ARRAY_BUFFER, vboPos);
            pglBufferData(GL_ARRAY_BUFFER, mesh.pos.size() * sizeof(float), mesh.pos.data(),
                          GL_STATIC_DRAW);

            pglGenBuffers(1, &vboCol);
            pglBindBuffer(GL_ARRAY_BUFFER, vboCol);
            pglBufferData(GL_ARRAY_BUFFER, mesh.nvert * 3, nullptr, GL_DYNAMIC_DRAW);

            pglGenBuffers(1, &ibo);
            pglBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
            pglBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh.tri.size() * sizeof(unsigned),
                          mesh.tri.data(), GL_STATIC_DRAW);

            pglBindBuffer(GL_ARRAY_BUFFER, 0);
            pglBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        }
    }
    void uploadColors() {
        if (!mesh.ok() || colors.empty())
            return;
        if (vboCol && pglBufferSubData) {
            pglBindBuffer(GL_ARRAY_BUFFER, vboCol);
            pglBufferSubData(GL_ARRAY_BUFFER, 0, colors.size(), colors.data());
            pglBindBuffer(GL_ARRAY_BUFFER, 0);
        }
    }

    bool load(int idx) {
        if (idx < 0 || idx >= (int)files.size())
            return false;
        if (idx == 0 || mesh.nvert == 0) {
            Mesh m;
            if (!parseVTK_full(files[idx], m))
                return false;
            mesh = std::move(m);
            globalMin = mesh.smin;
            globalMax = mesh.smax;
            setupGLBuffersOnce();
        } else {
            float smin, smax;
            if (!parseVTK_scalarsOnly(files[idx], mesh.nvert, mesh.scal, smin, smax))
                return false;
            mesh.smin = smin;
            mesh.smax = smax;
            if (!autoPerFrame) {
                globalMin = std::min(globalMin, smin);
                globalMax = std::max(globalMax, smax);
            }
        }
        rebuildColors();
        uploadColors();
        cur = idx;
        return true;
    }

    void rebuildColors() {
        if (!mesh.ok())
            return;
        colors.resize(mesh.nvert * 3);
        float smin = autoPerFrame ? mesh.smin : globalMin;
        float smax = autoPerFrame ? mesh.smax : globalMax;
        if (smax <= smin)
            smax = smin + 1.f;
        for (int i = 0; i < mesh.nvert; ++i) {
            float t = (mesh.scal[i] - smin) / (smax - smin);
            float r, g, b;
            if (cmap == GRAY)
                mapGray(t, r, g, b);
            else {
                float s = 2.f * t - 1.f;
                float t2 = 0.5f * (s + 1.f);
                mapRdBu(t2, r, g, b);
            }
            colors[i * 3 + 0] = (unsigned char)std::clamp((int)(r * 255), 0, 255);
            colors[i * 3 + 1] = (unsigned char)std::clamp((int)(g * 255), 0, 255);
            colors[i * 3 + 2] = (unsigned char)std::clamp((int)(b * 255), 0, 255);
        }
    }

    void next(int step) {
        int n = (int)files.size();
        if (n <= 1)
            return;
        int i = (cur + step) % n;
        if (i < 0)
            i += n;
        load(i);
    }
    void reset() {
        yaw = 0.f;
        pitch = 0.f;
        dist = 2.5f;
    }

    void draw() {
        glViewport(0, 0, winW, winH);
        glClearColor(bgWhite ? 1.f : 0.f, bgWhite ? 1.f : 0.f, bgWhite ? 1.f : 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        double aspect = winW > 0 ? double(winW) / double(std::max(1, winH)) : 1.0;
        double fovy = 45.0, zn = 0.1, zf = 100.0;
        double top = zn * std::tan(fovy * 0.5 * 3.141592653589793 / 180.0), right = top * aspect;
        glFrustum(-right, right, -top, top, zn, zf);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glTranslatef(0, 0, -dist);
        glRotatef(pitch * 57.29578f, 1, 0, 0);
        glRotatef(yaw * 57.29578f, 0, 1, 0);

        if (wire) {
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glDisable(GL_CULL_FACE);
        } else {
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glEnable(GL_CULL_FACE);
            glCullFace(GL_BACK);
        }

        if (!mesh.ok())
            return;

        if (hasVBO()) {
            pglBindBuffer(GL_ARRAY_BUFFER, vboPos);
            glEnableClientState(GL_VERTEX_ARRAY);
            glVertexPointer(3, GL_FLOAT, 0, (void *)0);

            pglBindBuffer(GL_ARRAY_BUFFER, vboCol);
            glEnableClientState(GL_COLOR_ARRAY);
            glColorPointer(3, GL_UNSIGNED_BYTE, 0, (void *)0);

            pglBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
            glDrawElements(GL_TRIANGLES, (GLsizei)mesh.tri.size(), GL_UNSIGNED_INT, (void *)0);

            glDisableClientState(GL_COLOR_ARRAY);
            glDisableClientState(GL_VERTEX_ARRAY);
            pglBindBuffer(GL_ARRAY_BUFFER, 0);
            pglBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        } else {
            glEnableClientState(GL_VERTEX_ARRAY);
            glVertexPointer(3, GL_FLOAT, 0, mesh.pos.data());
            glEnableClientState(GL_COLOR_ARRAY);
            glColorPointer(3, GL_UNSIGNED_BYTE, 0, colors.data());
            glDrawElements(GL_TRIANGLES, (GLsizei)mesh.tri.size(), GL_UNSIGNED_INT,
                           mesh.tri.data());
            glDisableClientState(GL_COLOR_ARRAY);
            glDisableClientState(GL_VERTEX_ARRAY);
        }
    }
};

// ---------- file dialog (1つだけ定義) ----------
static bool openFileDialog(std::string &outPath) {
#ifdef _WIN32
    char file[MAX_PATH] = {0};
    OPENFILENAMEA ofn{};
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = NULL;
    ofn.lpstrFile = file;
    ofn.nMaxFile = MAX_PATH;
    ofn.lpstrFilter = "VTK Files\0*.vtk\0All Files\0*.*\0";
    ofn.nFilterIndex = 1;
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST | OFN_EXPLORER;
    if (GetOpenFileNameA(&ofn)) {
        outPath = file;
        return true;
    }
#endif
    return false;
}

static void usage() {
    std::cout
        << "Usage:\n"
           "  vtk_viewer --file path\\to\\surface.vtk [--yawdeg 30 --pitchdeg -15 --dist 2.5]\n"
           "  vtk_viewer --series \"path\\prefix_%04d.vtk\" --first 0 --last 119 [--fps 24]\n"
           "  vtk_viewer --auto path\\to\\surface_0123.vtk  (自動で連番検出)\n"
           "  (引数なし→Windowsはファイルダイアログ／.vtkドロップ可／Oキーで開く)\n";
}

// ---------- main ----------
int main(int argc, char **argv) {
    std::string file, series;
    bool useAuto = false;
    int first = 0, lastFrame = -1;
    float yawdeg = 0, pitchdeg = 0, dist = 2.5f, fps = 24.f;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--file" && i + 1 < argc)
            file = argv[++i];
        else if (a == "--series" && i + 1 < argc)
            series = argv[++i];
        else if (a == "--first" && i + 1 < argc)
            first = std::atoi(argv[++i]);
        else if (a == "--last" && i + 1 < argc)
            lastFrame = std::atoi(argv[++i]);
        else if (a == "--auto" && i + 1 < argc) {
            file = argv[++i];
            useAuto = true;
        } else if (a == "--yawdeg" && i + 1 < argc)
            yawdeg = (float)std::atof(argv[++i]);
        else if (a == "--pitchdeg" && i + 1 < argc)
            pitchdeg = (float)std::atof(argv[++i]);
        else if (a == "--dist" && i + 1 < argc)
            dist = (float)std::atof(argv[++i]);
        else if (a == "--fps" && i + 1 < argc)
            fps = (float)std::atof(argv[++i]);
        else if (a == "-h" || a == "--help") {
            usage();
            return 0;
        }
    }

#ifdef _WIN32
    if (file.empty() && series.empty()) {
        std::string chosen;
        if (openFileDialog(chosen))
            file = chosen;
        else
            usage();
    }
#endif

    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_EVENTS) != 0) {
        std::cerr << "SDL_Init: " << SDL_GetError() << "\n";
        return 1;
    }
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    SDL_Window *win =
        SDL_CreateWindow("VTK Viewer (Fast)", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1200,
                         900, SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
    SDL_GLContext glc = SDL_GL_CreateContext(win);
    SDL_GL_SetSwapInterval(1); // 重いなら 0 に
    loadVBOFuncs();

    Viewer V;
    V.fps = fps;
    V.yaw = yawdeg * 0.017453292f;
    V.pitch = pitchdeg * 0.017453292f;
    V.dist = dist;

    if (!series.empty()) {
        for (int k = first; k <= lastFrame; ++k) {
            char buf[1024];
            std::snprintf(buf, sizeof(buf), series.c_str(), k);
            V.files.emplace_back(buf);
        }
        if (!V.load(0)) {
            std::cerr << "Failed to load first frame\n";
            return 1;
        }
        V.playing = (V.files.size() > 1);
    } else if (!file.empty()) {
        if (useAuto) {
            V.buildAutoSeries(file);
            V.playing = (V.files.size() > 1);
        } else
            V.files = {file};
        if (!V.load(0))
            std::cerr << "Failed to load " << V.files[0] << "\n";
    }

    bool running = true, dragging = false;
    int lastX = 0, lastY = 0;
    uint32_t lastTicks = SDL_GetTicks();
    double acc = 0.0;

    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT)
                running = false;
            else if (e.type == SDL_WINDOWEVENT && e.window.event == SDL_WINDOWEVENT_SIZE_CHANGED) {
                V.winW = e.window.data1;
                V.winH = e.window.data2;
            } else if (e.type == SDL_MOUSEBUTTONDOWN && e.button.button == SDL_BUTTON_LEFT) {
                dragging = true;
                lastX = e.button.x;
                lastY = e.button.y;
            } else if (e.type == SDL_MOUSEBUTTONUP && e.button.button == SDL_BUTTON_LEFT) {
                dragging = false;
            } else if (e.type == SDL_MOUSEMOTION && dragging) {
                int dx = e.motion.x - lastX, dy = e.motion.y - lastY;
                lastX = e.motion.x;
                lastY = e.motion.y;
                V.yaw += dx * 0.005f;
                V.pitch = std::clamp(V.pitch + dy * 0.005f, -1.5f, 1.5f);
            } else if (e.type == SDL_MOUSEWHEEL) {
                if (e.wheel.y > 0)
                    V.dist = std::max(0.5f, V.dist * 0.9f);
                if (e.wheel.y < 0)
                    V.dist = std::min(10.0f, V.dist * 1.1f);
            } else if (e.type == SDL_DROPFILE) {
                std::string p = e.drop.file ? e.drop.file : "";
                SDL_free(e.drop.file);
                if (!p.empty()) {
                    V.buildAutoSeries(p);
                    if (V.load(0))
                        V.playing = (V.files.size() > 1);
                }
            } else if (e.type == SDL_KEYDOWN) {
                switch (e.key.keysym.sym) {
                case SDLK_ESCAPE:
                    running = false;
                    break;
                case SDLK_SPACE:
                    V.playing = !V.playing;
                    break;
                case SDLK_LEFT:
                    V.next(-1);
                    break;
                case SDLK_RIGHT:
                    V.next(+1);
                    break;
                case SDLK_c:
                    V.cmap = (V.cmap == Viewer::GRAY ? Viewer::RDBU : Viewer::GRAY);
                    V.rebuildColors();
                    V.uploadColors();
                    break;
                case SDLK_a:
                    V.autoPerFrame = !V.autoPerFrame;
                    V.rebuildColors();
                    V.uploadColors();
                    break;
                case SDLK_b:
                    V.bgWhite = !V.bgWhite;
                    break;
                case SDLK_w:
                    V.wire = !V.wire;
                    break;
                case SDLK_r:
                    V.reset();
                    break;
                case SDLK_EQUALS:
                    V.fps = std::min(120.f, V.fps * 1.25f);
                    break;
                case SDLK_MINUS:
                    V.fps = std::max(1.f, V.fps * 0.8f);
                    break;
                case SDLK_o: {
                    std::string path;
                    if (openFileDialog(path)) {
                        V.buildAutoSeries(path);
                        if (V.load(0))
                            V.playing = (V.files.size() > 1);
                    }
                    break;
                }
                }
            }
        }
        uint32_t now = SDL_GetTicks();
        uint32_t dt = now - lastTicks;
        lastTicks = now;
        if (V.playing && V.files.size() > 1) {
            acc += dt / 1000.0;
            double dur = 1.0 / std::max(1.0f, V.fps);
            while (acc >= dur) {
                V.next(+1);
                acc -= dur;
            }
        }
        V.draw();
        SDL_GL_SwapWindow(win);
    }

    if (V.vboPos)
        pglDeleteBuffers(1, &V.vboPos);
    if (V.vboCol)
        pglDeleteBuffers(1, &V.vboCol);
    if (V.ibo)
        pglDeleteBuffers(1, &V.ibo);
    SDL_GL_DeleteContext(glc);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}
