#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#define N 60000
#define PI 3.1415926535   //圆周率
#define R 6378137         //地球半径
#define LAT 111319.55     //一纬度相差111319米
#define SGM 74.3          //标准差 64.7
#define BETA 34.56991141  // 0.06991141
using namespace std;

class myhash {
public:
    std::size_t operator()(const std::pair<int, int>& p) const {
        return std::hash<int>()(p.first) ^ std::hash<int>()(p.second);
    }
};

struct segment {
    double a, b, c;                     //线段方程 ax + by + c = 0
    pair<double, double> dot_a, dot_b;  //线段端点
    double length;
};

struct edge {
    /**
     **道路结构体
     **存储一个道路所有信息
     **/
    int id;                //路段id
    int p1, p2;            //路段起点id和终点id
    string way_string;     //路段等级名称
    int way_type;          //路段等级
    int c;                 //路段由c个地理采样点连起来构成
    vector<double> pot_x;  // c个地理采样点的纬度
    vector<double> pot_y;  // c个地理采样点的经度
    double length;         //路段长度
    vector<segment> line;  //一个路段是由几条线段连接而成
};

int m, n;
vector<vector<int>> vertice(N);                //与每个顶点相连的路段，邻接表存路网
vector<edge> my_edge;                          //所有的路段
vector<pair<double, double>> track;            //记录一条轨迹的所有采样点
vector<int> part[2500][2500];                  //将路网分块
vector<vector<pair<double, double>>> candDot;  //记录每个GPS采样点对应的候选点
vector<vector<int>> candDotline;               //记录每个候选点位于线段的第几段
vector<vector<int>> candEdgId;                 //记录每个候选路段的id
vector<vector<double>> matrixB;                //观测概率矩阵
vector<vector<double>> matrixA;                //状态转移矩阵
vector<vector<double>> matrixC;                //记录当前位置所有候选点的概率
vector<vector<int>> pre;
vector<int> ans;
unordered_map<pair<int, int>, double, myhash> path;  //记录两个顶点之间最短路
double dist[N];
int vis[N];

inline void clear_all() {
    //将上一条轨迹的全部状态信息清空
    candDot.clear();
    candDotline.clear();
    candEdgId.clear();
    matrixA.clear();
    matrixB.clear();
    matrixC.clear();
    track.clear();
    pre.clear();
    ans.clear();
}

inline void alloc_memory(unsigned int size_need) {
    /*******************
     * 给数组分配内存
     * 参数是需要的数组大小
     ********************/
    candDot.resize(size_need + 5);
    candDotline.resize(size_need + 5);
    candEdgId.resize(size_need + 5);
    matrixB.resize(size_need + 5);
    matrixC.resize(size_need + 5);
    pre.resize(size_need + 5);
    ans.resize(size_need + 5);
}

inline double dot_distance(const pair<double, double>& p1, const pair<double, double>& p2) {
    //计算两点之间的距离
    return sqrt((p1.first - p2.first) * (p1.first - p2.first) +
        (p1.second - p2.second) * (p1.second - p2.second));
}

inline segment get_seg(const pair<double, double>& dot_a, const pair<double, double>& dot_b) {
    double a = dot_a.second - dot_b.second;
    double b = dot_b.first - dot_a.first;
    double c = dot_a.first * dot_b.second - dot_b.first * dot_a.second;
    return { a, b, c, dot_a, dot_b, dot_distance(dot_a, dot_b) };
}

inline void add_part(const segment& seg, int id) {
    //将路段加入分块中
    int row_a = (seg.dot_a.first - 30.6) / 0.0005, col_a = (seg.dot_a.second - 121 * cos(PI / 6)) / 0.0005;
    int row_b = (seg.dot_b.first - 30.6) / 0.0005, col_b = (seg.dot_b.second - 121 * cos(PI / 6)) / 0.0005;
    int r1 = min(row_a, row_b), r2 = max(row_a, row_b), c1 = min(col_a, col_b), c2 = max(col_a, col_b);
    int ok = (seg.a * seg.b > 0);
    for (int r = r1; r <= r2; ++r) {
        for (int c = c1; c <= c2; ++c) {
            if (!ok) {
                double x1 = (r + 1) * 0.0005 + 30.6, y1 = c * 0.0005 + 121 * cos(PI / 6);
                double x2 = r * 0.0005 + 30.6, y2 = (c + 1) * 0.0005 + 121 * cos(PI / 6);
                if ((seg.a * x1 + seg.b * y1 + seg.c) * (seg.a * x2 + seg.b * y2 + seg.c) < 0) {
                    int in = 0;
                    for (int i = 0; i < part[r][c].size(); ++i) {
                        if (part[r][c][i] == id) {
                            in = 1;
                        }
                    }
                    if (!in)
                        part[r][c].push_back(id);
                }
            }
            else {
                double x1 = r * 0.0005 + 30.6, y1 = c * 0.0005 + 121 * cos(PI / 6);
                double x2 = (r + 1) * 0.0005 + 30.6, y2 = (c + 1) * 0.0005 + 121 * cos(PI / 6);
                if ((seg.a * x1 + seg.b * y1 + seg.c) * (seg.a * x2 + seg.b * y2 + seg.c) < 0) {
                    int in = 0;
                    for (int i = 0; i < part[r][c].size(); ++i) {
                        if (part[r][c][i] == id) {
                            in = 1;
                        }
                    }
                    if (!in)
                        part[r][c].push_back(id);
                }
            }
        }
    }
}

inline void read_edge() {
    int id, p1, p2, way_type, c;
    double x, y, length = 0;
    string way_string;
    vector<double> pot_x, pot_y;
    vector<segment> line;
    scanf("%d%d%d", &id, &p1, &p2);
    cin >> way_string;
    scanf("%d%d", &way_type, &c);
    for (int i = 1; i <= c; ++i) {
        scanf("%lf%lf", &x, &y);
        pot_x.push_back(x);
        pot_y.push_back(y * cos(PI / 6));
    }
    for (int i = 1; i < c; ++i) {
        line.push_back(get_seg({ pot_x[i - 1], pot_y[i - 1] }, { pot_x[i], pot_y[i] }));
    }
    for (int i = 0; i < line.size(); ++i) {
        length += line[i].length;
    }
    my_edge.push_back({ id, p1, p2, way_string, way_type, c, pot_x, pot_y, length, line });
    vertice[p1].push_back(id);
    for (int i = 0; i < line.size(); ++i) {
        add_part(line[i], id);
    }
}

inline void read_track() {
    /**************
     * 读入一条轨迹
     ***************/
    clear_all();
    int tim;
    double x, y;
    while ((scanf("%d", &tim), tim) > 10000) {
        scanf("%lf%lf", &x, &y);
        track.push_back({ x, y * cos(PI / 6) });
    }
    alloc_memory(track.size());
}

inline double shortest_path(int x, int y) {
    //找到x,y之间的最短路
    if (path.find({ x, y }) == path.end()) {
        priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> q;
        for (int i = 1; i < N; ++i) {
            dist[i] = 1e9;
            vis[i] = 0;
        }
        dist[x] = 0;
        q.push({ 0, x });
        while (!q.empty() && !vis[y]) {
            double d = q.top().first;
            int t = q.top().second;
            q.pop();
            if (vis[t])
                continue;
            dist[t] = d;
            vis[t] = 1;
            if (d > 0.05)
                continue;
            for (int i = 0; i < vertice[t].size(); ++i) {
                if (d + my_edge[vertice[t][i]].length < dist[my_edge[vertice[t][i]].p2]) {
                    dist[my_edge[vertice[t][i]].p2] = d + my_edge[vertice[t][i]].length;
                    q.push({ dist[my_edge[vertice[t][i]].p2], my_edge[vertice[t][i]].p2 });
                }
            }
        }
        path[{ x, y }] = dist[y];
    }
    return path[{ x, y }];
}

inline double dist_to_lend(int order, int x) {
    //候选点到路段左端点的距离
    double d = 0;
    int i;
    for (i = 0; i + 1 < candDotline[order][x]; ++i) {
        d += my_edge[candEdgId[order][x]].line[i].length;
    }
    d += dot_distance(my_edge[candEdgId[order][x]].line[i].dot_a, candDot[order][x]);
    return d;
}

inline double dist_to_rend(int order, int x) {
    //候选点到路段右端点的距离
    return my_edge[candEdgId[order][x]].length - dist_to_lend(order, x);
}

inline double real_distance(int order, int i, int j) {
    // 得到两个真实点之间的距离
    if (candEdgId[order - 1][i] == candEdgId[order][j]) {
        return abs(dist_to_lend(order, j) - dist_to_lend(order - 1, i));
    }
    return dist_to_rend(order - 1, i) + dist_to_lend(order, j) +
        shortest_path(my_edge[candEdgId[order - 1][i]].p2, my_edge[candEdgId[order][j]].p1);
}

inline pair<double, double> get_vector(const pair<double, double>& p1, const pair<double, double>& p2) {
    //通过两点计算出向量
    return { p2.first - p1.first, p2.second - p1.second };
}

inline double inner_product(const pair<double, double>& v1, const pair<double, double>& v2) {
    // 计算向量内积
    return v1.first * v2.first + v1.second * v2.second;
}

inline pair<double, double> get_pedal(const pair<double, double>& dot, const segment& seg) {
    //找到点到线段上的垂足
    pair<double, double> vec_ap = get_vector(seg.dot_a, dot);
    pair<double, double> vec_bp = get_vector(seg.dot_b, dot);
    pair<double, double> vec_ab = get_vector(seg.dot_a, seg.dot_b);
    double r = inner_product(vec_ap, vec_ab) / pow(dot_distance(seg.dot_a, seg.dot_b), 2);
    if (r <= 0)
        return seg.dot_a;
    if (r >= 1)
        return seg.dot_b;
    double x0 = dot.first, y0 = dot.second;
    double a = seg.a, b = seg.b, c = seg.c;
    double x = (b * b * x0 - a * b * y0 - a * c) / (a * a + b * b);
    double y = (a * a * y0 - a * b * x0 - b * c) / (a * a + b * b);
    return { x, y };
}

inline pair<double, double> get_canddot(int order, const edge& edg) {
    //找到候选点并将候选点在路段的第几段上记录下来
    pair<double, double> ans = { 9999.0, 9999.0 }, t;
    int num = 0;
    for (int i = 0; i < edg.line.size(); ++i) {
        t = get_pedal(track[order], edg.line[i]);
        if (dot_distance(track[order], t) < dot_distance(track[order], ans)) {
            ans = t;
            num = i;
        }
    }
    candDotline[order].push_back(num);
    return ans;
}

inline void find_cand(int order) {
    //找到所有候选点
    int cnt[130000] = { 0 };
    int row = (track[order].first - 30.6) / 0.0005, col = (track[order].second - 121 * cos(PI / 6)) / 0.0005;
    for (int i = max(row - 1, 0); i <= min(row + 1, 2499); ++i) {
        for (int j = max(col - 1, 0); j <= min(col + 1, 2499); ++j) {
            for (int k = 0; k < part[i][j].size(); ++k) {
                int pos = part[i][j][k];
                if (cnt[pos] == 0) {
                    pair<double, double> tmp = get_canddot(order, my_edge[pos]);
                    // if (dot_distance(tmp, track[order]) * LAT < 70.0) {
                    candEdgId[order].push_back(pos);
                    candDot[order].push_back(tmp);
                    //}
                    cnt[pos] = 1;
                }
            }
        }
    }
    if (candDot[order].size() == 0) {
        candDot[order] = candDot[order - 1];
        candEdgId[order] = candEdgId[order - 1];
        candDotline[order] = candDotline[order - 1];
    }
}

inline void get_matA(int order) {
    /******************************************
     * 得到状态转移矩阵
     * 参数是一个下标，
     * 表示得到第x-1个状态到第x个状态的状态转移矩阵
     ********************************************/
    int row = candDot[order - 1].size(), col = candDot[order].size();
    matrixA.clear();
    matrixA.resize(row);
    vector<double> beta;
    for (int i = 0; i < row; ++i) {
        matrixA[i].resize(col);
        for (int j = 0; j < col; ++j) {
            double d = LAT * abs(dot_distance(track[order - 1], track[order]) - real_distance(order, i, j));
            matrixA[i][j] = exp(-d / BETA / BETA);
            // beta.push_back(d);
            // matrixA[i][j] = d;
        }
    }
    // sort(beta.begin(), beta.end());
    // double b;
    // if (beta.size() % 2) b = beta[beta.size() / 2] / log(2);
    // else b = (beta[beta.size() / 2] + beta[beta.size() / 2 - 1]) / 2 / log(2);
    // for (int i = 0; i < row; ++i) {
    //     for (int j = 0; j < col; ++j) {
    //         matrixA[i][j] = exp(-matrixA[i][j] / b) / b;
    //     }
    // }
}

inline void get_matB(int order) {
    find_cand(order);
    vector<double> sigma;
    double s = 0;
    for (int i = 0; i < candDot[order].size(); ++i) {
        double d = LAT * dot_distance(candDot[order][i], track[order]);
        double probability = exp(-d * d / (2 * SGM * SGM)) / (sqrt(2 * PI) * SGM);
        matrixB[order].push_back(probability);
        // sigma.push_back(d);
        // matrixB[order].push_back(d);
    }
    // sort(sigma.begin(), sigma.end());
    // if (sigma.size() % 2) s = 1.4826 * sigma[sigma.size() / 2];
    // else s = 1.4826 * (sigma[sigma.size() / 2] + sigma[sigma.size() / 2 - 1]);
    // //printf("%lf", s);
    // for (int i = 0; i < matrixB[order].size(); ++i) {
    //     double d = matrixB[order][i];
    //     matrixB[order][i] = exp(-d * d / (2 * s * s)) / (sqrt(2 * PI) * s);
    // }
}

inline void get_matC(int order) {
    get_matB(order);
    if (order == 0) {
        matrixC[order] = matrixB[order];
        return;
    }
    get_matA(order);
    matrixC[order].resize(matrixB[order].size(), 0.0);
    pre[order].resize(matrixC[order].size(), 0);
    double t = 0;
    for (int i = 0; i < matrixB[order].size(); ++i) {
        for (int j = 0; j < matrixB[order - 1].size(); ++j) {
            if (matrixC[order][i] < matrixC[order - 1][j] + matrixB[order][i] * matrixA[j][i]) {
                matrixC[order][i] = matrixC[order - 1][j] + matrixB[order][i] * matrixA[j][i];
                pre[order][i] = j;
            }
            t = max(t, matrixC[order][i]);
        }
    }
    for (int i = 0; i < matrixC[order].size(); ++i) {
        matrixC[order][i] /= t;
    }
}

inline void find_road() {
    for (int i = 0; i < track.size(); ++i) {
        get_matC(i);
    }
    int end = 0;
    for (int i = 0; i < matrixC[track.size() - 1].size(); ++i) {
        if (matrixC[track.size() - 1][i] > matrixC[track.size() - 1][end]) {
            end = i;
        }
    }
    for (int i = track.size() - 1; i > 0; --i) {
        ans[i - 1] = pre[i][end];
        end = ans[i - 1];
    }
    for (int i = 0; i < track.size(); ++i) {
        printf("%d ", candEdgId[i][ans[i]]);
    }
    printf("\n");
}

int main() {
    // freopen("D:\\Code\\Code_c\\Single\\map_matching\\sample.in", "r", stdin);
    // freopen("D:\\Code\\Code_c\\Single\\map_matching\\1.out", "w", stdout);
    scanf("%d", &n);
    for (int i = 1; i <= n; ++i) {
        read_edge();
    }
    scanf("%d", &m);
    printf("%d\n", m);
    for (int i = 1; i <= m; ++i) {
        read_track();
        find_road();
    }
    return 0;
}