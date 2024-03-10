#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <filesystem> // C++17 or later

#define MAX_GRID 10000
#define precision 0.0005
#define STEP 0.1
#define STEP_angle 0.01
#define epsilon precision
#define origin Point{0,0}
#define input_localminmax "C:\\Users\\sebas\\OneDrive\\Desktop\\EPFL\\Odd area\\local_minimums.txt"
#define ROTATE 0
#define output_localminmax (ROTATE == 0 ? "C:\\Users\\sebas\\OneDrive\\Desktop\\EPFL\\Odd area\\oddarea_localminmax_shift_up.txt": "C:\\Users\\sebas\\OneDrive\\Desktop\\EPFL\\Odd area\\oddarea_localminmax_rotate.txt")


struct Point {
    double x;
    double y;
    // Overload the + operator for Point structs
    Point operator+(const Point& other) const {
        return {x + other.x, y + other.y};
    }
};

Point rotate (Point p, double alpha){
    //rotates clockwise by alpha
    return {p.x * cos(alpha) + p.y * sin(alpha), -p.x * sin(alpha) + p.y * cos(alpha)};
}

// Overloading the << operator for Point
std::ostream& operator<<(std::ostream& os, const Point& point) {
    os << "(" << point.x << ", " << point.y << ")";
    return os;
}


// Function to calculate the angle of a point with respect to the reference point
double calculateAngle(const Point& reference, const Point& point) {
    return atan2(point.x - reference.x, point.y - reference.y);
}

// Custom comparison function for sorting points in anticlockwise order
// starting at the point {-1,0}
bool comparePoints(const Point& reference, const Point& a, const Point& b) {
    return calculateAngle(reference, a) < calculateAngle(reference, b);
}



class OddArea {
public:
    OddArea(std::vector<Point> centers) : centers(centers), nr_circles(0),
                                                                   grid_top(-MAX_GRID), grid_bottom(MAX_GRID), grid_left(MAX_GRID), grid_right(-MAX_GRID) {
        nr_circles = centers.size();
        for (const auto &center : centers) {
            grid_top = std::max(grid_top, center.y + 1);
            grid_bottom = std::min(grid_bottom, center.y - 1);
            grid_left = std::min(grid_left, center.x - 1);
            grid_right = std::max(grid_right, center.x + 1);
        }

    }

    double dist(Point point1, Point point2) {
        return std::sqrt(std::pow(point1.x - point2.x, 2) + std::pow(point1.y - point2.y, 2));
    }

    int is_odd(const Point& point) {
        int ans = 0;
        for (const auto &center : centers) {
            ans += dist(point, center) < 1;
        }
        return ans % 2;
    }

    double calculate_odd_area() {
        int odd_points = 0;
        for (double i = grid_left; i < grid_right; i += precision) {
            for (double j = grid_bottom; j < grid_top; j += precision) {
                Point point = {i, j};
                odd_points += is_odd(point);
            }
        }
        this -> odd_area = odd_points * std::pow(precision, 2);
        return this -> odd_area;
    }
    Point calculate_gradient(int center_index){
        //calculates gradient with the sum from the paper by Pinchasi
        Point center = centers[center_index];
        std::vector<Point> intersections;
        for (int i=0; i<nr_circles; ++i){
            double dist_center = dist(center, centers[i]);
            if (i == center_index ||  dist_center>=2 ) continue;

            double x1 = center.x, y1 = center.y;
            double x2 = centers[i].x, y2 = centers[i].y;

            Point intersection1 = {(x1+x2)/2 + (y2-y1)/2*std::sqrt(4/(std::pow(dist_center,2))-1),
                                    (y1+y2)/2 +(x1-x2)/2*std::sqrt(4/(std::pow(dist_center,2))-1)};
            Point intersection2 = {(x1+x2)/2 - (y2-y1)/2*std::sqrt(4/(std::pow(dist_center,2))-1),
                                    (y1+y2)/2 - (x1-x2)/2*std::sqrt(4/(std::pow(dist_center,2))-1)};;

            intersections.push_back(intersection1);
            intersections.push_back(intersection2);
        }

        // Sort the points in anticlockwise order with respect to the center
        std::sort(intersections.begin(), intersections.end(), [&center](const Point& a, const Point& b) {
            return comparePoints(center, a, b);
        });

//      Uncomment to print the intersections
//        for (const auto &intersection : intersections){
//            std::cout<<"intersection: "<<intersection<<"atan: "<< atan2(intersection.x,intersection.y)<<std::endl;
//        }

        Point sum = {0,0};
        int flip =-1;
        for(const auto &intersection: intersections){
            sum.x += 2*flip*intersection.x;
            sum.y += 2*flip*intersection.y;
            flip*= -1;
        }

        return {sum.y, -sum.x};
    }

    Point calculate_gradient_direct(int center_index){
        // calculates gradient with the definition by computing the odd area of a small shift
        std::vector<Point> new_centers_x(centers.size()), new_centers_y(centers.size());
        std::copy(centers.begin(), centers.end(), new_centers_x.begin());
        std::copy(centers.begin(), centers.end(), new_centers_y.begin());
        new_centers_x[center_index].x += epsilon;
        new_centers_y[center_index].y += epsilon;
        double current_odd_area = this->calculate_odd_area();
        OddArea new_circles_x(new_centers_x), new_circles_y(new_centers_y);

        return {double((new_circles_x.calculate_odd_area()- current_odd_area)/epsilon),
                double((new_circles_y.calculate_odd_area()- current_odd_area)/epsilon)};
    }


    void one_step_GD(int center_index){
        Point gradient = this->calculate_gradient(center_index);
        this->centers[center_index] = {this->centers[center_index].x -gradient.x*STEP, this->centers[center_index].y -gradient.y*STEP};
        std::cout<< "previous odd area "<<this->odd_area<<std::endl;
        this->calculate_odd_area();
        std::cout<< "odd area now"<<this->odd_area<<std::endl;
    }

    void all_step_GD(){
        // This function moves along the gradient all ceneters one by one in the general setting
        Point gradient;
        for (int i=0; i< this->nr_circles; ++i){
            gradient = this->calculate_gradient(i);
            this->centers[i] = {this->centers[i].x -gradient.x*STEP, this->centers[i].y -gradient.y*STEP};
        }
        std::cout<< "previous odd area "<<this->odd_area<<std::endl;
        this->calculate_odd_area();
        std::cout<< "odd area now "<<this->odd_area<<std::endl;
    }

    void check_local_min (){
        Point unit = {0, 0.005};
        Point initial = this->centers[0];
        double initial_odd_area = this->calculate_odd_area();
        std::cout <<"Odd area:"<< initial_odd_area << std::endl;
        for (int i=0; i<10; ++i){
            this->centers[0] = initial + rotate(unit, 2*M_PI*i/10);
            double curr_odd_area = this -> calculate_odd_area();
            std::cout <<"Move by "<<rotate(unit, 2*M_PI*i/10) <<"-> odd area initial minus current:"<< initial_odd_area - curr_odd_area << std::endl;
        }
    }

    void savePoints(const std::string& filename) {
            std::ofstream outputFile(filename);

            if (outputFile.is_open()) {
                for (const auto center: centers) {
                    outputFile << center.x << "," << center.y << "\n";
                }

                outputFile.close();
                std::cout << "Points saved to " << filename << std::endl;
            } else {
                std::cerr << "Unable to open file: " << filename << std::endl;
            }
        }

    std::vector<Point> centers;
    int nr_circles;
    double grid_top;
    double grid_bottom;
    double grid_left;
    double grid_right;
    double odd_area = 0;

};

class OddArea_concentric :public OddArea{
public:

    OddArea_concentric(std::vector<Point> centers): OddArea (centers)
        {radius = dist(origin, centers[0]);};
    //to be initialized with points centered around the origin {0,0}

    void one_step_GD(int center_index){
        // now we have to stay on the circle when we move along the gradient
        Point center = this->centers[center_index];
        Point gradient = this->calculate_gradient(center_index);
        Point tangent = {-center.y, center.x};

        double new_angle = calculateAngle(origin, center);

        if (tangent.x*gradient.x + tangent.y*gradient.y < 0){
            //anticlockwise
             new_angle -= STEP_angle;
        }
        else{
            //clockwise
           new_angle += STEP_angle;
        }
        this->centers[center_index] = {this->radius*cos(new_angle), this->radius*sin(new_angle) };
        this->centers[center_index] = {this->centers[center_index].x -gradient.x*STEP, this->centers[center_index].y -gradient.y*STEP};
        std::cout<< "previous odd area "<<this->odd_area<<std::endl;
        this->calculate_odd_area();
        std::cout<< "odd area now "<<this->odd_area<<std::endl;

    }
    double radius;
};

std::vector<Point> generateEquidistantPoints(int n, double radius) {
    std::vector<Point> points;

    for (int i = 0; i < n; ++i) {
        double theta = M_PI/2 + 2.0 * M_PI * i / n;
        double x = radius * cos(theta);
        double y = radius * sin(theta);

        points.push_back({x, y});
    }

    return points;
}

void check_local_minmax (){
    std::ifstream file(input_localminmax);
    double n, r;
    std::ofstream outputFile(output_localminmax);
    if (ROTATE)
        outputFile << "Rotating clockwise by angle 0.01 \n\n";
    else
        outputFile << "Shifting by adding {0.005, 0.005} \n\n";
    int count = 0, total =0;
    while (file >> n >> r) {
        total++;
        std::vector<Point> centers = generateEquidistantPoints(n, r);
        OddArea circles(centers);
        outputFile << "n = "<< n <<", r = "<<r<<std::endl;
        outputFile<<circles.centers[0]<<" has gradient "<<circles.calculate_gradient(0)<<std::endl;
        double prev_area = circles.calculate_odd_area();
        outputFile << "Odd area: "<< prev_area << std::endl;
        double curr_area = 0;
        if (ROTATE == 1){
            circles.centers[0] = rotate(circles.centers[0], STEP_angle);
            curr_area = circles.calculate_odd_area();
            outputFile << "Odd area after rotation: "<< curr_area << std::endl;}
        else{
            circles.centers[0] = circles.centers[0] + Point({0, 0.005});
            curr_area = circles.calculate_odd_area();
            outputFile << "Odd area after shift: "<< curr_area << std::endl;
        }
        if (curr_area > prev_area){
            count +=1;
            outputFile<< "STOP, the diff is "<<curr_area -prev_area<<std::endl<<std::endl;
        }
        outputFile << std::endl;
        }

        outputFile << "We have found " <<count <<"out of "<< total;
}


int main() {

    std::vector<Point> centers = generateEquidistantPoints(13, 1.0632816909765221);
   // std::vector<Point> centers = {{0,0},{1,0},{0.5,sqrt(3)/2}};
    OddArea circles(centers);
   // circles.savePoints("C:\\Users\\sebas\\OneDrive\\Desktop\\EPFL\\Odd area\\Centers_before.csv");
  //  std::cout << circles.calculate_odd_area() << std::endl;
   // std::cout<<circles.centers[0]<<" has gradient "<<circles.calculate_gradient(0)<<std::endl;
  //  std::cout<<circles.centers[0]<<" has gradient direct "<<circles.calculate_gradient_direct(0)<<std::endl;

   // circles.centers[0] = rotate(circles.centers[0], STEP_angle);
   // std::cout << circles.calculate_odd_area() << std::endl;
  //,  circles.check_local_min();
//    std::cout<<circles.calculate_gradient_direct(0) <<std::endl;

//    for (int i=0;i<200;++i)
//        circles.one_step_GD(i%10);
    circles.savePoints("C:\\Users\\sebas\\OneDrive\\Desktop\\EPFL\\Odd area\\Centers_after.csv");
 // check_local_minmax();
    return 0;
}
