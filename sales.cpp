#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <random>

struct Coord {
  double lon = 0.0;
  double lat = 0.0;
};

struct Params {
  // outputs
  std::string inFile;
  std::string outRoute = "route.dat";
  std::string outSched = "anneal.csv";


  // SA controls
  int steps = 3000000;
  double T0 = 5000.0;
  double alpha = 0.999997;
  int printEvery = 200000;
  int schedEvery = 2000;

  // move mix
  double swapProb = 0.05;

  // reproducibility
  uint64_t seed = 0; // 0 => use time()
};

static void Die(const std::string& msg) {
  std::cerr << msg << "\n";
  std::exit(1);
}

static Params ParseArgs(int argc, char** argv) {
  Params p;
  if (argc < 2) {
    Die("Usage: tsp_sa_alt cities.dat [routeout.dat] [schedule.csv] "
        "[--steps N] [--T0 X] [--alpha X] [--printEvery N] [--schedEvery N] "
        "[--swapProb X] [--seed N]");
  }

  // positional args
  p.inFile = argv[1];
  if (argc > 2 && argv[2][0] != '-') p.outRoute = argv[2];
  if (argc > 3 && argv[3][0] != '-') p.outSched = argv[3];

  // optional flags
  for (int i = 2; i < argc; ++i) {
    std::string a = argv[i];
    auto need = [&](const char* name) -> std::string {
      if (i + 1 >= argc) Die(std::string("Missing value for ") + name);
      return std::string(argv[++i]);
    };

    if (a == "--steps") p.steps = std::stoi(need("--steps"));
    else if (a == "--T0") p.T0 = std::stod(need("--T0"));
    else if (a == "--alpha") p.alpha = std::stod(need("--alpha"));
    else if (a == "--printEvery") p.printEvery = std::stoi(need("--printEvery"));
    else if (a == "--schedEvery") p.schedEvery = std::stoi(need("--schedEvery"));
    else if (a == "--swapProb") p.swapProb = std::stod(need("--swapProb"));
    else if (a == "--seed") p.seed = static_cast<uint64_t>(std::stoull(need("--seed")));
  }

  if (p.steps <= 0) Die("Error: --steps must be positive.");
  if (p.T0 <= 0.0) Die("Error: --T0 must be positive.");
  if (p.alpha <= 0.0 || p.alpha >= 1.0) Die("Error: --alpha must be in (0,1).");
  if (p.swapProb < 0.0 || p.swapProb > 1.0) Die("Error: --swapProb must be in [0,1].");
  if (p.printEvery <= 0) p.printEvery = 200000;
  if (p.schedEvery <= 0) p.schedEvery = 2000;

  return p;
}

static std::vector<Coord> ReadCities(const std::string& fname, int NMAX = 3000) {
  std::ifstream in(fname);
  if (!in) Die("Error: cannot open " + fname);

  std::vector<Coord> cities;
  cities.reserve(1024);

  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    if (static_cast<int>(cities.size()) >= NMAX) break;

    std::istringstream iss(line);
    Coord c;
    if (iss >> c.lon >> c.lat) cities.push_back(c);
  }

  if (cities.size() < 2) Die("Error: need at least 2 cities.");
  return cities;
}

static double HaversineKm(double lat1, double lon1, double lat2, double lon2) {
  const double R = 6371.0;
  const double rad = M_PI / 180.0;

  const double la1 = lat1 * rad;
  const double la2 = lat2 * rad;
  const double dlat = (lat2 - lat1) * rad;
  const double dlon = (lon2 - lon1) * rad;

  const double s1 = std::sin(dlat * 0.5);
  const double s2 = std::sin(dlon * 0.5);
  const double a = s1 * s1 + std::cos(la1) * std::cos(la2) * (s2 * s2);
  const double c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
  return R * c;
}

static double TourLengthKm(const std::vector<Coord>& route) {
  const int n = static_cast<int>(route.size());
  double sum = 0.0;
  for (int i = 0; i < n - 1; ++i) {
    sum += HaversineKm(route[i].lat, route[i].lon, route[i + 1].lat, route[i + 1].lon);
  }
  sum += HaversineKm(route[n - 1].lat, route[n - 1].lon, route[0].lat, route[0].lon);
  return sum;
}

class Annealer {
 public:
  explicit Annealer(uint64_t seed)
      : rng_(seed),
        u01_(0.0, 1.0) {}

  void SwapTwo(std::vector<Coord>& r) {
    int n = (int)r.size();
    std::uniform_int_distribution<int> pick(0, n - 1);
    int i = pick(rng_);
    int j = pick(rng_);
    while (j == i) j = pick(rng_);
    std::swap(r[i], r[j]);
  }

  void ReverseSegment(std::vector<Coord>& r) {
    int n = (int)r.size();
    // choose [a,b] inclusive
    std::uniform_int_distribution<int> pickA(0, n - 2);
    int a = pickA(rng_);
    std::uniform_int_distribution<int> pickB(a + 1, n - 1);
    int b = pickB(rng_);
    std::reverse(r.begin() + a, r.begin() + b + 1);
  }

  bool Accept(double dE, double T) {
    if (dE <= 0.0) return true;
    double p = std::exp(-dE / T);
    return (u01_(rng_) < p);
  }

  double U01() { return u01_(rng_); }

 private:
  std::mt19937_64 rng_;
  std::uniform_real_distribution<double> u01_;
};

static void WriteRoute(const std::string& fname, const std::vector<Coord>& best) {
  std::ofstream out(fname);
  if (!out) Die("Error: cannot write route file " + fname);

  out << std::fixed << std::setprecision(6);
  for (const auto& c : best) out << c.lon << " " << c.lat << "\n";
}

int main(int argc, char** argv) {
  Params P = ParseArgs(argc, argv);

  uint64_t seed = P.seed;
  if (seed == 0) seed = static_cast<uint64_t>(std::time(nullptr));

  std::vector<Coord> cur = ReadCities(P.inFile, 3000);
  std::vector<Coord> trial = cur;
  std::vector<Coord> best = cur;

  std::cout << "Read " << cur.size() << " cities from data file\n";
  std::cout << std::fixed;

  double Ecur = TourLengthKm(cur);
  double Ebest = Ecur;

  std::cout << "Initial distance (km): " << std::setprecision(3) << Ecur << "\n";

  std::ofstream sched(P.outSched);
  if (!sched) Die("Error: cannot write schedule file " + P.outSched);
  sched << "T,current_km,best_km\n";

  Annealer A(seed);

  double T = P.T0;
  auto t0 = std::chrono::steady_clock::now();

  for (int s = 0; s < P.steps; ++s) {
    trial = cur;

    // Choose move
    if (A.U01() < P.swapProb) A.SwapTwo(trial);
    else A.ReverseSegment(trial);

    double Etrial = TourLengthKm(trial);
    double dE = Etrial - Ecur;

    if (A.Accept(dE, T)) {
      cur.swap(trial);
      Ecur = Etrial;
      if (Ecur < Ebest) {
        Ebest = Ecur;
        best = cur;
      }
    }

    T *= P.alpha;

    if (s % P.schedEvery == 0) {
      sched << std::setprecision(12) << T << ","
            << std::setprecision(12) << Ecur << ","
            << std::setprecision(12) << Ebest << "\n";
    }

    if (s % P.printEvery == 0) {
      std::cout << "step " << s
                << "  T=" << std::setprecision(3) << T
                << "  best=" << std::setprecision(3) << Ebest << "\n";
    }
  }

  auto t1 = std::chrono::steady_clock::now();
  double sec = std::chrono::duration<double>(t1 - t0).count();

  WriteRoute(P.outRoute, best);

  std::cout << "Best distance (km):    " << std::setprecision(3) << Ebest << "\n";
  std::cout << "Runtime (s):           " << std::setprecision(3) << sec << "\n";
  std::cout << "Wrote route:           " << P.outRoute << "\n";
  std::cout << "Wrote schedule:        " << P.outSched << "\n";

  return 0;
}

