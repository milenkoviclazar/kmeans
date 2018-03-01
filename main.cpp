#include <iostream>
#include <vector>
#include <random>
#include <cassert>

using namespace std;

const double EPSILON = 1.0;

random_device rd;
mt19937 gen(rd());

struct Point {
    double x, y;

    Point() {
        x = 0.0;
        y = 0.0;
    }

    Point(double _x, double _y) : x(_x), y(_y) {}
};

double dist(const Point &a, const Point &b) {
    return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}

vector<int> findNearestCenters(const vector<Point> &coords, const vector<Point> &centers) {
    vector<int> nearest(coords.size(), 0);
    for (int i = 0; i < coords.size(); i++) {
        nearest[i] = 0;
        for (int j = 1; j < centers.size(); j++) {
            if (dist(coords[i], centers[j]) < dist(coords[i], centers[nearest[i]])) {
                nearest[i] = j;
            }
        }
    }
    return nearest;
}

double calculateCost(const vector<Point> &coords, const vector<Point> centers) {
    vector<int> nearest = findNearestCenters(coords, centers);
    double cost = 0.0;
    for (int i = 0; i < coords.size(); i++) {
        cost += dist(coords[i], centers[nearest[i]]);
    }
    return cost;
}

double calculateCost(const vector<Point> &coords, const vector<Point> centers, const vector<int> nearest) {
    double cost = 0.0;
    for (int i = 0; i < coords.size(); i++) {
        cost += dist(coords[i], centers[nearest[i]]);
    }
    return cost;
}

void lloyd(const vector<Point> &coords, vector<Point> &centers) {
    int n = coords.size();
    int k = centers.size();
    vector<int> nearest = findNearestCenters(coords, centers);
    double oldCost = calculateCost(coords, centers, nearest);
    double newCost;
    vector<int> counter(centers.size(), 0);
    do {
        fill(centers.begin(), centers.end(), Point(0, 0));
        fill(counter.begin(), counter.end(), 0);
        for (int i = 0; i < coords.size(); i++) {
            centers[nearest[i]].x += coords[i].x;
            centers[nearest[i]].y += coords[i].y;
            counter[nearest[i]]++;
        }
        for (int i = 0; i < centers.size(); i++) {
            centers[i].x /= counter[i];
            centers[i].y /= counter[i];
        }
        nearest = findNearestCenters(coords, centers);
        newCost = calculateCost(coords, centers, nearest);
    } while (fabs(oldCost - newCost) > EPSILON);
}

void kMeansPlusPlus(int k, const vector<Point> &coords) {
    vector<Point> centers;
    int n = coords.size();
    uniform_int_distribution<> die(0, coords.size() - 1);
    centers.push_back(coords[die(gen)]);
    vector<double> weights(coords.size(), 0.0);
    for (int i = 0; i < coords.size(); i++) {
        weights[i] = dist(coords[i], centers[0]);
    }
    int nCandidates = max(1, (int)log(k));
    vector<int> candidates(nCandidates, 0);
    for (int i = 1; i < k; i++) {
        // D^2 sampling
        discrete_distribution<int> weightedDistribution(weights.begin(), weights.end());
        for (int j = 0; j < nCandidates; j++) {
            candidates[j] = weightedDistribution(gen);
        }
        int bestIdx = -1;
        double bestGain = -1.0;
        for (int j = 0; j < nCandidates; j++) {
            centers.push_back(coords[candidates[j]]);
            double gain = calculateCost(coords, centers);
            if (bestIdx == -1 || gain < bestGain) {
                bestGain = gain;
                bestIdx = j;
            }
            centers.pop_back();
        }
        assert(bestIdx != -1);
        centers.push_back(centers[candidates[bestIdx]]);
        // update the weights
        for (int j = 0; j < coords.size(); i++) {
            weights[j] = min(weights[j], dist(coords[j], centers[i]));
        }

    }
    lloyd(coords, centers);
}

void kMeansGreedy(int k, const vector<Point> &coords) {
    vector<Point> centers;
    int n = coords.size();
    for (int i = 0; i < k; i++) {
        int bestIdx = -1;
        double bestGain = -1.0;
        for (int j = 0; j < n; j++) {
            // TODO: check if unselected
            centers.push_back(coords[j]);
            double gain = calculateCost(coords, centers);
            if (bestIdx == -1 || gain < bestGain) {
                bestGain = gain;
                bestIdx = j;
            }
            centers.pop_back();
        }
        assert(bestIdx != -1);
        centers.push_back(coords[bestIdx]);
    }
    lloyd(coords, centers);
}

int main(int argc, char *argv[]) {
    int n, k;
    vector<Point> coords;
    if (argc == 1) {
        cerr << "Invalid number of arguments. Set the first argument to be:\n"
                "'random' for random generation of points\n"
                "'file' for reading from the file" << endl;
        return 0;
    }
    if (string(argv[1]) == "random") {
        uniform_int_distribution<> nRandom(10, 1000);
        n = nRandom(gen);
        uniform_int_distribution<> kRandom(1, n);
        k = kRandom(gen);
        uniform_int_distribution<> coordRandom(-10000, 10000);
        for (int i = 0; i < n; i++) {
            coords.push_back(Point(coordRandom(gen), coordRandom(gen)));
        }
    } else if (string(argv[1]) == "file") {
        cin >> n >> k;
        for (int i = 0; i < n; i++) {
            int x, y;
            cin >> x >> y;
            coords.push_back(Point(x, y));
        }
    } else {
        cerr << "Invalid argument" << endl;
        return 0;
    }
    kMeansPlusPlus(k, coords);
    kMeansGreedy(k, coords);
    return 0;
}