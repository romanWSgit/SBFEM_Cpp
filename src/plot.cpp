//
// Created by Roman Wallner- Silberhuber on 04.06.23.
//

#include "plot.h"
#include "helper_functions.h"

void plotSbfemSuperelement (const SuperElement &sE, bool gnuplot)
{
    if (gnuplot)
        {

//            Gnuplot gp ("/opt/local/bin/gnuplot"); // This is a shitty idea !!! Path is hardcoded
            Gnuplot gp;
            // set the terminal to Aqua
            gp << "set terminal aqua\n";

            Eigen::MatrixXd Data = sE.getMSEElementsM ();

            Eigen::MatrixXd mat_extended (Data.rows (), Data.cols () + 2);
            mat_extended
                << Data, Eigen::MatrixXd::Constant (Data.rows (), 1, 0.0), Eigen::MatrixXd::Constant (Data.rows (), 1,
                                                                                                      0.0);
            Data = mat_extended;



            // Color of triangles
            std::string color = "0xff0000";
            gp << "set title 'Superelement'  font 'Helvetica Bold,20'\n";
            // Enable grid with custom line style and color
            gp << "set grid lc rgb 'gray' lt 1 lw 1\n";
            gp << "set linetype 1 lc rgb \"black\"\n";
            for (int i = 0; i < Data.rows (); ++i)
                {
                    // creating the polygon for each triangle in the data
                    gp << "set object " << i + 1 << " polygon from "
                       << Data (i, 0) << "," << Data (i, 1) << " to "
                       << Data (i, 2) << "," << Data (i, 3) << " to "
                       << Data (i, 4) << "," << Data (i, 5) << " to "
                       << Data (i, 0) << "," << Data (i, 1) << "\n";

                    // setting the color and border for each polygon
                    gp << "set object " << i + 1 << " fc rgb \"" << color
                       << "\" fs solid 0.5 border lt 1 lw 3\n";

                    // calculate triangle centroid
                    double centroid_x =
                        (Data (i, 0) + Data (i, 2) + Data (i, 4)) / 3;
                    double centroid_y =
                        (Data (i, 1) + Data (i, 3) + Data (i, 5)) / 3;

                    // Set the label for each polygon
                    gp << "set label " << i + 1 << " \"" << i + 1 << "\" at "
                       << centroid_x << "," << centroid_y \
 << "center font \"Helvetica Bold,16\"\n";
                }

            gp << "set xrange[-13:13]\n";
            gp << "set yrange[-6:6]\n";
            gp << "set size ratio -1\n";
            gp << "plot NaN notitle\n";
            gp << "pause mouse close\n";
        }
    else
        {
            Eigen::MatrixXd elementsM = sE.getMSEElementsM ();
            const int iElem = elementsM.rows ();

            std::vector<double> centre = {0.0, 0.0}; // scaling centre
            std::cout << centre[0] << " " << centre[1] << std::endl;

//        std::array<std::array<double, 4>, iElem> elements = {{
//            {{-12.0, -3.0, 0.0, -3.0}},
//            {{0.0, -3.0, 12.0, -3.0}}
//        }};


            matplot::hold (matplot::on);
            for (int i = 0; i < iElem; ++i)
                {
                    std::vector<double> x = {centre[0], elementsM (i, 0),
                                             elementsM (i, 2)};
                    std::vector<double> y = {centre[1], elementsM (i, 1),
                                             elementsM (i, 3)};
                    double centroidX = (x[0] + x[1] + x[2]) / 3;
                    double centroidY = (y[0] + y[1] + y[2]) / 3;
                    x.push_back (x[0]); // Schließen des Dreiecks
                    y.push_back (y[0]);

                    matplot::fill (x, y, "r");

                    matplot::plot (x, y, "b")->line_width (3);

                    matplot::text (centroidX, centroidY,
                                   "el_{" + std::to_string (i + 1) + "}")
                        ->font_size (16);

                }
            std::vector<double> centreVx = {
                centre[0]}; // X-Koordinate des Punktes
            std::vector<double> centreVy = {
                centre[1]}; // Y-Koordinate des Punktes
            //matplot::plot(centreVx,centreVy, "o")->line_width(3);// 'o' ist der Markerstil für einen runden Punkt
//   matplot::text(centre[0], centre[1], "centre");
            matplot::textarrow (1, 1, centre[0], centre[1], "centre");
            matplot::title ("Superelement");
            matplot::axis (matplot::equal);
            matplot::grid (matplot::on);
            matplot::show ();
        }
}

void visualizeDoubleMatrix (const Eigen::MatrixXd &matrix)
{
    std::vector<std::vector<double>> data =
        eigenToStdVector2D (matrix);
    {
        using namespace matplot;

        auto h = figure(true);
        h->name("Diagonal Visualization");
        h->number_title(false);
        h->color("green");
        h->position (350,30);
        h->size(1000, 1000);
        h->draw();
        h->font("Arial");
        h->font_size(15);
        h->title("Diagonal of T");
        heatmap(data);
        colormap(palette::rdbu ());
        auto ax1 = gca();
        ax1->axes_aspect_ratio(1.0);
        h->draw();

        show();
    }

}

void visualizeOnlyDiagonalOfDoubleMatrix (const Eigen::MatrixXd &matrix)
{
    Eigen::MatrixXd diagonalM = matrix.diagonal().asDiagonal();

    std::vector<std::vector<double>> data =
        eigenToStdVector2D (diagonalM);

    {
        using namespace matplot;

        auto h = figure(true);
        h->name("Diagonal Visualization");
        h->number_title(false);
        h->color("green");
        h->position (350,30);
        h->size(1000, 1000);
        h->draw();
        h->font("Arial");
        h->font_size(15);
        h->title("Diagonal of T");
        heatmap(data);
        colormap(palette::rdbu ());
        auto ax1 = gca();
        ax1->axes_aspect_ratio(1.0);
        h->draw();

        show();
    }
}