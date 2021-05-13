using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace kalman
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void button1_Click(object sender, EventArgs e)
        {

            Random r = new Random();

            List<Matrix<double>> F = new List<Matrix<double>>();
            List<Matrix<double>> B = new List<Matrix<double>>();
            List<Matrix<double>> P = new List<Matrix<double>>();
            List<Matrix<double>> H = new List<Matrix<double>>();
            List<Matrix<double>> R = new List<Matrix<double>>();
            List<Matrix<double>> K = new List<Matrix<double>>();

            List<Vector<double>> x = new List<Vector<double>>();
            List<Vector<double>> u = new List<Vector<double>>();
            List<Vector<double>> z = new List<Vector<double>>();

            List<Vector<double>> extrapolated_x = new List<Vector<double>>();
            List<Vector<double>> cor_x = new List<Vector<double>>();

            List<Matrix<double>> extrapolated_P = new List<Matrix<double>>();

            Matrix<double> I = Matrix<double>.Build.DenseIdentity(3);

            int delta_t = 1;
            int v_max = 10;
            int omega_max = 10;
            int count = 100;
            int x_max = 10;
            int y_max = 10;
            int x_min = 0;
            int y_min = 0;
            double sigma_x = (x_max - x_min)/6;
            double sigma_y = (y_max - y_min) / 6;
            double sigma_teta = 1.097;
            int x0 = 0;
            int y0 = 0;
            int ksi = 10;
            int tay = 10;
            double percent = 0.1;
            double r_min = Math.Sqrt((ksi - x_min)*(ksi - x_min) + (tay - y_min)*(tay - y_min));
            double r_max = Math.Sqrt((ksi - x_max) * (ksi - x_max) + (tay - y_max) * (tay - y_max));
            double sigma_r = (r_max - r_min) / 6;
            double sigma_phi = 1.097;
            double[] teta = new double[count];
            double[] v = new double[count];
            double[] omega = new double[count];
            v[0] = r.NextDouble();
            omega[0] = r.NextDouble();
            for (int i = 1; i < count; i++)
            {
                v[i] = r.NextDouble() * 2 * v[i - 1] * percent - v[i - 1] * percent;         // TODO 
                if (v[i] > v_max)
                {
                    count--;
                    continue;
                }
                
            }
            for (int i = 0; i < count; i++)
            {
                teta[i] = r.NextDouble() * Math.PI * 2 - Math.PI; // radians
            }
            for (int i = 1; i < count; i++)
            {
                omega[i] = r.NextDouble() * 2 * omega[i - 1] * percent - omega[i - 1] * percent; // TODO
                if (Math.Abs(v[i]) > omega_max)
                {
                    count--;
                    continue;
                }
                Console.WriteLine(i + " " + omega[i]);
            }
            // наблюдаемый x
            x.Add(Vector<double>.Build.DenseOfArray(new double[] { x0, y0, teta[0] }));
            for (int i = 1; i < count; i++)
            {
                x.Add(Vector<double>.Build.DenseOfArray(new double[] 
                {   r.NextDouble() * x_max,
                    r.NextDouble() * y_max,
                    teta[i]
                }));
            }
            // u
            for (int i = 0; i < count; i++)
            {
                u.Add(Vector<double>.Build.DenseOfArray(new double[]
                {   v[i],
                    omega[i]
                }));
            }
            // z
            for (int i = 0; i < count; i++)
            {
                z.Add(Vector<double>.Build.DenseOfArray(new double[]
                {   Math.Sqrt((ksi - x[i][0])*(ksi - x[i][0]) + (tay - x[i][1])*(tay - x[i][1])),
                    Math.Atan2((tay - x[i][1]), (ksi - x[i][0])) - teta[i]
                }));
            }
            // F
            for (int i = 0; i < count; i++)
            {
                F.Add(Matrix<double>.Build.DenseOfArray(new double[,]
                {
                    { 1, 0, -v[i] * Math.Sin(teta[i] * delta_t) },
                    { 0, 1, v[i] * Math.Cos(teta[i] * delta_t) },
                    { 0, 0, 1 }
                }));
            }

            // B
            for (int i = 0; i < count; i++)
            {
                B.Add(Matrix<double>.Build.DenseOfArray(new double[,]
                {
                    { Math.Cos(teta[i] * delta_t), 0},
                    { Math.Sin(teta[i] * delta_t), 0 },
                    { 0, delta_t }
                }));
            }

            // P
            //for (int i = 0; i < count; i++)
            //{
            P.Add(Matrix<double>.Build.DenseOfArray(new double[,]
            {
                { sigma_x*sigma_x, 0, 0 },
                { 0, sigma_y*sigma_y, 0 },
                { 0, 0, sigma_teta }
            }));
            //}

            // H
            for (int i = 0; i < count; i++)
            {
                H.Add(Matrix<double>.Build.DenseOfArray(new double[,]
                {
                    { -2*(ksi - x[i][0])/Math.Sqrt((ksi - x[i][0]) * (ksi - x[i][0]) + (tay - x[i][1]) * (tay - x[i][1])),
                      -2*(tay - x[i][1])/Math.Sqrt((ksi - x[i][0]) * (ksi - x[i][0]) + (tay - x[i][1]) * (tay - x[i][1])),
                       0 },
                    {  (tay - x[i][1])/((ksi - x[i][0]) * (ksi - x[i][0]) + (tay - x[i][1]) * (tay - x[i][1])),
                       -1*(ksi - x[i][0]) * (ksi - x[i][0])/((ksi - x[i][0]) * (ksi - x[i][0]) + (tay - x[i][1]) * (tay - x[i][1])),
                      -1 }
                }));
            }

            // R
            for (int i = 0; i < count; i++)
            {
                R.Add(Matrix<double>.Build.DenseOfArray(new double[,]
                {
                    { sigma_r*sigma_r, 0 },
                    { 0, sigma_phi }
                }));
            }
            extrapolated_x.Add(x[0]);
            extrapolated_P.Add(P[0]);
            K.Add(Matrix<double>.Build.DenseOfArray(new double[,]
                {
                    { 0, 0 },
                    { 0, 0 }
                }));
            cor_x.Add(x[0]);
            for (int i = 1; i < count; i++)
            {
                // 1. Экстраполяция состояния
                extrapolated_x.Add(
                    F[i] * x[i - 1] + B[i] * u[i]
                    );
                // 2. Экстраполяция матрицы ковариации
                extrapolated_P.Add(
                    F[i] * P[i - 1] * F[i].Transpose() //+ P[0]
                    );
                // 3. Усиление ко Калману
                K.Add(
                    extrapolated_P[i] * H[i].Transpose() * (H[i] * extrapolated_P[i] * H[i].Transpose() + R[i]).Inverse()
                    );
                // 4. Коррекция вектора состояния
                cor_x.Add(
                    extrapolated_x[i] + K[i] * (z[i] - H[i] * extrapolated_x[i]) 
                    );
                // 5.Расчёт матрицы ковариации
                P.Add(
                    (I - K[i] * H[i]) * extrapolated_P[i]
                    );
            }
            //List<Vector<double>> output = new List<Vector<double>>();
            DataTable output = new DataTable("Kalman_result");
            // создаем столбцы для таблицы Kalman_result

            DataColumn ex_x      = new DataColumn("X", Type.GetType("System.Double"));
            DataColumn ex_y      = new DataColumn("Y", Type.GetType("System.Double"));
            DataColumn ex_teta   = new DataColumn("Teta", Type.GetType("System.Double"));
            DataColumn corr_x    = new DataColumn("X^", Type.GetType("System.Double"));
            DataColumn corr_y    = new DataColumn("Y^", Type.GetType("System.Double"));
            DataColumn corr_teta = new DataColumn("Teta^", Type.GetType("System.Double"));
           // DataColumn ob_x      = new DataColumn("^X", Type.GetType("System.Double"));
           // DataColumn ob_y      = new DataColumn("^Y", Type.GetType("System.Double"));
           // DataColumn ob_teta   = new DataColumn("^Teta", Type.GetType("System.Double"));

            //output.Columns.Add(ob_x);
            output.Columns.Add(ex_x);
            output.Columns.Add(corr_x);
            //output.Columns.Add(ob_y);
            output.Columns.Add(ex_y);
            output.Columns.Add(corr_y);
            //output.Columns.Add(ob_teta);
            output.Columns.Add(ex_teta);
            output.Columns.Add(corr_teta);

            DataRow row = output.NewRow();
            for (int i = 0; i < count; i++)
            {
                output.Rows.Add(new object[]
                {
                   // x[i][0],              // ^X
                    extrapolated_x[i][0], // X
                    cor_x[i][0],          // X^
                    //x[i][1],              // ^Y
                    extrapolated_x[i][1], // Y
                    cor_x[i][1],          // Y^
                    //x[i][2],              // ^Teta
                    extrapolated_x[i][2], // Teta
                    cor_x[i][2]           // Teta^
                }); // добавляем строку
            }

            dataGridView1.DataSource = output;
        }
    }
}






/*Console.WriteLine($"M {M}");
           Console.WriteLine($"V {V}");
           Console.WriteLine($"M*V {MV}");
           Console.WriteLine($"V*M {VM}");*/

/*

Matrix<double> F = Matrix<double>.Build.DenseOfArray(new double[,]
       {
                { 1, 2 },
                { 3, 6 }
       });

Matrix<double> B = Matrix<double>.Build.DenseOfArray(new double[,]
{
                { 1, 2 },
                { 3, 6 }
});

Vector<double> x = Vector<double>.Build.DenseOfArray(new double[] { 3, 4 });
Vector<double> u = Vector<double>.Build.DenseOfArray(new double[] { 3, 4 });

    */


/*
       StringBuilder stb1 = new StringBuilder();
       stb1.Append("[");
       foreach (var pt in h)
       {
           stb1.Append("x:" + pt[0] + "y:" + pt[1] + " , ");
       }
       stb1.Append("]");
       string result1 = stb1.ToString();*/
//textBox2.Text = result;