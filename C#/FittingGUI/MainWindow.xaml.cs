using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Controls.Primitives;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

using ParFilt;
using DirectoryNotFoundException = System.IO.DirectoryNotFoundException;


namespace FittingGUI
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        
        public MainWindow()
        {
            CultureInfo ci = new CultureInfo("en-EN");
            Thread.CurrentThread.CurrentCulture = ci;
            Thread.CurrentThread.CurrentUICulture = ci;


            Fs = 44100;
            maxPoles = 48;
            FixedNPoles = 48;
            NFIR = 1;
            MyPlot = new Plot();


            InitializeComponent();
            DataContext = this;
        }

        private void DragCompleted(object sender, DragCompletedEventArgs e)
        {
            ComputeResponse();
        }

        private void SliderChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            var Slider = sender as Slider;
            if (Slider.Template != null) //slider exists
            {
                var track = Slider.Template.FindName("PART_Track", Slider) as Track;
                if (track == null || track.Thumb.IsDragging) //slider has not null track
                {
                    if (FixedNPolesCheck.IsChecked ?? false) //fixed n is checked -> change npolesHigh slider value
                    {
                        NPolesHighSlider.Value = FixedNPoles - NPolesLowSlider.Value;
                    }

                    return;
                }
            }
            else
            {
                return;
            }

            ComputeResponse();
        }

        private void FixedNpolesChange(object sender, RoutedEventArgs e)
        {
            if (FixedNPolesCheck.IsChecked ?? false)
            {
                MaxPolesBox.IsEnabled = false;
                NPolesLowSlider.Maximum = FixedNPoles;
                NPolesHighSlider.Maximum = FixedNPoles;
                NPolesHighSlider.Value = FixedNPoles - NPolesLowSlider.Value;
                NPolesHighSlider.IsEnabled = false;
                NPolesHighBox.IsEnabled = false;
            }
            else
            {
                MaxPolesBox.IsEnabled = true;
                NPolesLowSlider.Maximum = MaxPoles;
                NPolesHighSlider.Maximum = MaxPoles;
                NPolesHighSlider.IsEnabled = true;
                NPolesHighBox.IsEnabled = true;
            }

        }

        private void FixedPolesBoxTextChanged(object sender, TextChangedEventArgs e)
        {
            NPolesLowSlider.Maximum = FixedNPoles;
            NPolesHighSlider.Maximum = FixedNPoles;
            NPolesLowSlider.Value = Math.Min(NPolesLowSlider.Value, FixedNPoles);
            NPolesHighSlider.Value = FixedNPoles - NPolesLowSlider.Value;
        }

        private void ComputeResponse()
        {
            var crossFr = CrossFreqSlider.Value;
            var crossLen = Convert.ToInt32(CrossLengthFreqSlider.Value);
            var npoles1 = Convert.ToInt32(NPolesLowSlider.Value);
            var npoles2 = Convert.ToInt32(NPolesHighSlider.Value);
            var lambda1 = LambdaLowSlider.Value;
            var lambda2 = LambdaHighSlider.Value;


            if (target != null)
            {
                Array.Resize(ref filter, target.Length);
                int size = target.Length;
                ParFiltDesign.computeResponse(filter, W, target, out size, npoles1, npoles2, crossFr, crossLen,
                    lambda1, lambda2, Fs,NFIR,useNAK);
                MyPlot.AddToPlotdB(fr, filter, "filter", 1);
            }
        }



        private void LoadResponse(object sender, RoutedEventArgs e)
        {
            int size = 1600;
            double[] w = new double[size];
            double[] H = new double[size];
            

            bool result = false;
            if (ItemSelected == null)
                return;
            do
            {
                string ff = FilePath + "\\" + ItemSelected;

                result = ParFiltDesign.loadResponse(w,H, out size, FilePath + "\\" + ItemSelected);

                if (!result)
                {
                    if (size < 0)   //error occured
                        return;
                    w = new double[size];
                    H = new double[size];
                }
            } while (!result);

            Array.Resize(ref w, size); //shortens the output to correct length
            Array.Resize(ref H, size); //shortens the output to correct length

            fr = w.Select(r => r * Fs/(2*Math.PI)).ToArray();
            MyPlot.AddToPlotdB(fr, H, "Target", 0);
            target = H;
            W = w;
        }


        // PROPERTIES
        public int FixedNPoles { get; set; }
        public int Fs { get; set; }
        public int NFIR { get; set; }
        public bool useNAK { get; set; }
        public string ItemSelected { get; set; }

        public int MaxPoles //{ get; set; }
        {
            get => maxPoles;
            set
            {
                maxPoles = Math.Max(value, 1);
                NPolesLowSlider.Maximum = MaxPoles;
                NPolesHighSlider.Maximum = MaxPoles;
                NPolesLowSlider.Value = Math.Min(MaxPoles, NPolesLowSlider.Value);
                NPolesHighSlider.Value = Math.Min(MaxPoles, NPolesHighSlider.Value);
            }
        }

        public Plot MyPlot { get; }

        public string FilePath { get; set; }
        public string ExportFileName { get; set; }


        // PRIVATE FIELDS
        private double[] target;
        private double[] W;
        private double[] fr;
        private double[] filter;
        //private bool isDragging = false;
        private int maxPoles;

        private void EnterPressed(object sender, KeyEventArgs e)
        {
            if (e.Key != System.Windows.Input.Key.Enter) return;

            try
            {
                DirectoryInfo dinfo = new DirectoryInfo(FilePath);
                FileInfo[] files = dinfo.GetFiles("*.csv*");
                foreach (FileInfo file in files)
                {
                    FolderContentsListBox.Items.Add(file.Name);
                
                }
            }
            catch (Exception exception)
            {
                return;
            }
        }

        private void UseNAKCheckClick(object sender, RoutedEventArgs e)
        {
            ComputeResponse();
        }

        private void ExportCoefficientsToCsv(object sender, RoutedEventArgs e)
        {
            var crossFr = CrossFreqSlider.Value;
            var crossLen = Convert.ToInt32(CrossLengthFreqSlider.Value);
            var npoles1 = Convert.ToInt32(NPolesLowSlider.Value);
            var npoles2 = Convert.ToInt32(NPolesHighSlider.Value);
            var lambda1 = LambdaLowSlider.Value;
            var lambda2 = LambdaHighSlider.Value;
            var currPath = FilePathBox.Text;

            if (target != null && ExportFileName != "")
            {
                int size = target.Length;
                ParFiltDesign.exportToCsv(W,target,out size,npoles1,npoles2,crossFr,crossLen,lambda1,lambda2,Fs,NFIR,useNAK,currPath + "\\" + ExportFileName);
            }
        }

        private void InfoBtnMessage(object sender, RoutedEventArgs e)
        {
            MessageBox.Show(
                            "First enter path to location of stored files with frequency responses, press Enter\n"+
                            "You should see available files in provided folder. To load selected file press Load and frequency response will be plotted.\n\n" +
                            "Supported file format is csv and file should contain 2 columns, first angular frequencies, second magnitude of frequency response (not in dB).\n\n" +
                            "Filter design is based on Balazs Bank dual warping method. Frequency response is devided into 2 parts, each frequency vector is warped by lambda parameter and poles of filter are calculated using Kallman method." +
                            "Then iterative Steiglitz-McBride method optimizes the poles positions. Filter is divided into parallel second-order sections and zeros of filters are calculated.\n\n" +
                            "Inside the algorithms spline interpolation and extrapolation is used. You can choose between linear extrapolation or not-a-knot type.\n\n" +
                            "You can export the coefficients of parallel second-order sections and FIR part of the filter by Export button. File will be located at same directory as path.\n\n" +
                            "For more information see:\nhttp://home.mit.bme.hu/~bank/parfilt/\n\n" +
                            "App by Matouš Vrbík\n" +
                            "matousvrbik@gmail.com"
                ,
                "Information",
                MessageBoxButton.OK,
                MessageBoxImage.Information);
        }
    }
}
