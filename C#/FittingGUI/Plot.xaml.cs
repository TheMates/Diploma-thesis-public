using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using LiveCharts;
using LiveCharts.Configurations;
using LiveCharts.Defaults;
using LiveCharts.Geared;
using LiveCharts.Wpf;

namespace FittingGUI
{
    public partial class Plot : UserControl
    {
        public Plot()
        {
            InitializeComponent();

            var mapper = Mappers.Xy<ObservablePoint>()
                            .X(point => Math.Log10(point.X))
                            .Y(point => point.Y);

            SeriesCollection = new SeriesCollection(mapper);


            FormatterX = value => Math.Pow(10, value).ToString("N4");

            DataContext = this;
        }

        public void AddToPlot(IEnumerable<double> x, IEnumerable<double> y, string title, int index)
        {
            //modifying the series collection will animate and update the chart
            ChartValues<ObservablePoint> nove = new ChartValues<ObservablePoint>();
            ObservablePoint [] points = new ObservablePoint[x.Count()];
            for (var i = 0; i < x.Count(); ++i)
            {
                points[i] = new ObservablePoint { X = x.ElementAt(i), Y = y.ElementAt(i) };
            }
            nove.AddRange(points);


            SeriesCollection.Add(new LineSeries
            {
                Title = title,
                Values = nove.AsGearedValues().WithQuality(Quality.Medium) ,
                LineSmoothness = 0, //0: straight lines, 1: really smooth lines
                //PointGeometry = Geometry.Parse("m 25 70.36218 20 -28 -20 22 -8 -6 z"),
                PointGeometry = DefaultGeometries.Circle,
                PointGeometrySize = 1,
                //PointForeground = Brushes.Gray,
                Fill = Brushes.Transparent
            });
            
        }

        public void AddToPlotdB(IEnumerable<double> x, IEnumerable<double> y, string title, int index)
        {
            //modifying the series collection will animate and update the chart
            ChartValues<ObservablePoint> nove = new ChartValues<ObservablePoint>();
            ObservablePoint[] points = new ObservablePoint[x.Count()];
            for (var i = 0; i < x.Count(); ++i)
            {
                points[i] = new ObservablePoint { X = x.ElementAt(i), Y = 20*Math.Log10(y.ElementAt(i)) };
            }
            nove.AddRange(points);
            if (SeriesCollection.Count < index + 1)
            {
                SeriesCollection.Add(new LineSeries());
            }

            SeriesCollection.CurrentSeriesIndex = index;
            if (index == 0 && SeriesCollection.Count > 1)
                SeriesCollection.RemoveAt(1);

                SeriesCollection[index] = (new LineSeries
            {
                Title = title,
                Values = nove.AsGearedValues().WithQuality(Quality.Medium),
                LineSmoothness = 0, //0: straight lines, 1: really smooth lines
                //PointGeometry = Geometry.Parse("m 25 70.36218 20 -28 -20 22 -8 -6 z"),
                PointGeometry = DefaultGeometries.Circle,
                PointGeometrySize = 1,
                Fill = Brushes.Transparent
            });
            

        }


        public SeriesCollection SeriesCollection { get; set; }
        public string[] Labels { get; set; }

        public Func<double, string> FormatterX { get; set; }
        public Func<double, string> YFormatter { get; set; }

        private void ResetZoom(object sender, MouseButtonEventArgs e)
        {
            XAxis.MinValue = double.NaN;
            XAxis.MaxValue = double.NaN;
            YAxis.MinValue = double.NaN;
            YAxis.MaxValue = double.NaN;
        }    
    }
}