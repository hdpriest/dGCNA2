package networkCalcPackage;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.TreeMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class Operations {

    public static float GINI(float[] array1, float[] array2) {
        float GINI_coeff;
        float Numerator = 0.0f;
        float Denominator = 0.0f;
        int[] sortRanks1 = getIndicesInOrder(array1);
        int[] sortRanks2 = getIndicesInOrder(array2);
        for (int i = 0; i < array1.length; i++) {
            float v2 = ((2 * (i + 1)) - array1.length - 1) * array1[sortRanks1[i]];
            float v1 = ((2 * (i + 1)) - array1.length - 1) * array1[sortRanks2[i]];
            Denominator += v2;
            Numerator += v1;
        }
        GINI_coeff = Numerator / Denominator;
        return GINI_coeff;
    }

    public static GCNMatrix copyNames(String[] names, GCNMatrix NetB) {
        /// if you always copy Rows->rows and columns, you can't go wrong
        NetB.setRowNames(names);
        NetB.setAlternateRowNames(names);
        NetB.setColumnNames(names);
        return NetB;
    }

    public static int[] getRanks(float[] array) {
        Map<Integer, Float> map = new HashMap<Integer, Float>(array.length);
        for (int i = 0; i < array.length; i++) {
            map.put(i, array[i]);
        }

        List<Entry<Integer, Float>> l
                = new ArrayList<Entry<Integer, Float>>(map.entrySet());

        Collections.sort(l, new Comparator<Entry<?, Float>>() {
            @Override
            public int compare(Entry<?, Float> e1, Entry<?, Float> e2) {
                return e2.getValue().compareTo(e1.getValue());
            }
        });

        int[] result = new int[array.length];
        String input = "";
        String ranks = "";
        for (int i = 0; i < result.length; i++) {
            result[l.get(i).getKey()]=i+1;
        }
        return result;
    }
    
    public static int[] getIndicesInOrder(float[] array) {
        Map<Integer, Float> map = new HashMap<Integer, Float>(array.length);
        for (int i = 0; i < array.length; i++) {
            map.put(i, array[i]);
        }

        List<Entry<Integer, Float>> l
                = new ArrayList<Entry<Integer, Float>>(map.entrySet());

        Collections.sort(l, new Comparator<Entry<?, Float>>() {
            @Override
            public int compare(Entry<?, Float> e1, Entry<?, Float> e2) {
                return e2.getValue().compareTo(e1.getValue());
            }
        });

        int[] result = new int[array.length];
        String input = "";
        String ranks = "";
        for (int i = 0; i < result.length; i++) {
            result[i] = l.get(i).getKey();
        }
//        System.err.println(ranks+"\n"+input+"\n");
        return result;
    }

    public static float[] compareNetworksViaTOM(GCNMatrix Net1, GCNMatrix Net2) {
        int D = Net1.getNumRows();
        GCNMatrix ReturnFrame = new GCNMatrix(D, D);
        ReturnFrame.setRowNames(Net1.getRowNames());
        ReturnFrame.setAlternateRowNames(Net2.getRowNames());
        float[] cTOMs = new float[D];
        //ReturnFrame = Operations.copyNames(Net1.getRowNames(), ReturnFrame);
        for (int i = 0; i < D; i++) {
            for (int j = 0; j < D; j++) {
                float T = 0;
                if (i == j) {
                    float product = 0;
                    float i_k = Net1.findK(i, i);
                    float j_k = Net2.findK(j, j);
                    for (int u = 0; u < D; u++) {
                        if ((u != i) && (u != j) && (Net1.testValue(i, u)) && (Net2.testValue(j, u))) {
                            float i_v = Net1.getValueByEntry(i, u);
                            float j_v = Net2.getValueByEntry(j, u);
                            float max = Math.max(i_v, j_v);
                            product += i_v * j_v / max;
                            /// if node is not connected to anything, all products are zero
                        }
                    }
                    float k_min = Math.min(i_k, j_k);
                    float DFIJ = 0;
                    T = (product + DFIJ) / (k_min + 1 - DFIJ); // if one node unconnected, = 0+0/0+1-0
                    cTOMs[i] = T;
                } else {
                }

                //ReturnFrame.setValueByEntry(T, i, j);
            }
        }
        return cTOMs;
    }
    
    public static float[] compareNetworksViaAverage(GCNMatrix Net1, GCNMatrix Net2) {
        int D = Net1.getNumRows();
        float[] cTOMs = new float[D];
        for (int i = 0; i < D; i++) {
            for (int j = 0; j < D; j++) {
                float T = 0;
                if (i == j) {
                    float sum = 0.0f;
                    float i_k = Net1.findK(i, i);
                    float j_k = Net2.findK(j, j);
                    float a_k = 0.0f;
                    for (int u = 0; u < D; u++) {
                        if ((u != i) && (u != j) && ((Net1.testValue(i, u)) || (Net2.testValue(j, u)))) {
                            float i_v = Net1.getValueByEntry(i, u);
                            float j_v = Net2.getValueByEntry(j, u);
                            float max = Math.max(i_v, j_v);
                            sum += max * Math.abs(j_v - i_v);
                            a_k += max;
                        }
                    }
                    float avg = sum / a_k; 
                    cTOMs[i] = avg;
                } else {
                }

                //ReturnFrame.setValueByEntry(T, i, j);
            }
        }
        return cTOMs;
    }
    public static float[] compareNetworksViaMWW(GCNMatrix Net1, GCNMatrix Net2) {
        int D = Net1.getNumRows();
        float[] cTOMs = new float[D];
        for (int i = 0; i < D; i++) {
            int j=i;
            cTOMs[i]=1.0f;
            ArrayList<Double> Net1Values = new ArrayList<Double>();
            ArrayList<Double> Net2Values = new ArrayList<Double>();
            float a_k = 0.0f;
            for (int u = 0; u < D; u++) {
                if ((u != i) && (u != j) && ((Net1.testValue(i, u)) && (Net2.testValue(j, u)))) {
                    float i_v = Net1.getValueByEntry(i, u);
                    float j_v = Net2.getValueByEntry(j, u);
                    if((i_v<0.05) && (j_v <0.05)) continue;
                    Net1Values.add(Double.valueOf(i_v));
                    Net2Values.add(Double.valueOf(j_v));
                    a_k += 1.0f;
                }
            }
            if(a_k <= 100.0f)continue;
            double[] net1values = new double[Net1Values.size()];
            double[] net2values = new double[Net2Values.size()];
            int X =0;
            for(Double value : Net1Values){
                net1values[X]= (double) value;
                X++;
            }
            X=0;
            for(Double value : Net2Values){
                net2values[X]= (double) value;
                X++;
            }
            MannWhitneyUTest MWW = new MannWhitneyUTest();
            double PofU = MWW.mannWhitneyUTest(net1values, net2values);
            cTOMs[i] = (float) PofU;
        }
        return cTOMs;
    }

    public static GCNMatrix calculateTOM(GCNMatrix Adjacency, int Threads) {
        int D = Adjacency.getNumRows();
        GCNMatrix ReturnMatrix = new GCNMatrix(D, D);
        ReturnMatrix.setRowNames(Adjacency.getRowNames());
        ReturnMatrix.setAlternateRowNames(Adjacency.getAlternateRowNames());
        ReturnMatrix.setColumnNames(Adjacency.getColumnNames());
        ExecutorService pool = Executors.newFixedThreadPool(Threads);
        ExecutorCompletionService<HashMap<String, float[]>> completionService = new ExecutorCompletionService<>(pool);
        List<Future<HashMap<String, float[]>>> taskList = new ArrayList<Future<HashMap<String, float[]>>>();
        ConcurrentLinkedQueue<String> queue = new ConcurrentLinkedQueue<String>();
        //System.err.println("Processing topological overlap using " + Threads + " threads.");
        for (int i = 0; i < D; i++) {
            String S = String.valueOf(i);
            queue.add(S);
        }
        int Number = D * D;
        System.err.println("Added "+ Number +" Tasks To Multithreaded Processing Engine");
        for (int i = 0; i < Threads; i++) {
            Callable<HashMap<String, float[]>> worker = new ConcurrentProcessing(Adjacency, queue, "tom");
            Future<HashMap<String, float[]>> submit = completionService.submit(worker);
            taskList.add(submit);

        }
        
        for (int t = 0; t < Threads; t++) {
            try {
                HashMap<String, float[]> hm = completionService.take().get();
                //System.err.println("obtained result for thread " + t);
                int r = 0;
                for (Map.Entry<String, float[]> entry : hm.entrySet()) {
                    String s = entry.getKey();
                    int i = Integer.valueOf(s);
                    int size = D - i;
                    float[] d = new float[size];
                    d = entry.getValue();
                    for (int j = 0; j < d.length; j++) {
                        int coord = j + i;
                        ReturnMatrix.setValueByEntry(d[j], i, coord);
                        r++;
                    }
                }
                //System.err.println("Processed "+ r + " records");
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
            //System.err.println("Thread " + t + " complete.");
        }
        //System.err.println("Done.");
        pool.shutdownNow();
        return ReturnMatrix;
    }

    private static void explode() {
        System.err.println("No correlation got to this point!?");
        System.exit(0);
    }
    
    public static GCNMatrix scaleMatrix(GCNMatrix OldMatrix, float scaleFactor, int Threads) {
        // prepwork
        int D = OldMatrix.getNumRows();
        GCNMatrix NewMatrix = new GCNMatrix(D, D);
        NewMatrix.setRowNames(OldMatrix.getRowNames());
        NewMatrix.setAlternateRowNames(OldMatrix.getAlternateRowNames());
        NewMatrix.setColumnNames(OldMatrix.getRowNames());
        // thread prep
        ExecutorService pool = Executors.newFixedThreadPool(Threads);
        ExecutorCompletionService<HashMap<String, float[]>> completionService = new ExecutorCompletionService<>(pool);
        List<Future<HashMap<String, float[]>>> taskList = new ArrayList<Future<HashMap<String, float[]>>>();
        ConcurrentLinkedQueue<String> queue = new ConcurrentLinkedQueue<String>();

        //queue prep
        for (int i = 0; i < D; i++) {
            String S = String.valueOf(i);
            queue.add(S);
        }
        int Number = D * D;
        System.err.println("Added "+ Number +" Tasks To Multithreaded Processing Engine");
        //add the tasks
        for (int i = 0; i < Threads; i++) {
            Callable<HashMap<String, float[]>> worker = new ConcurrentProcessing(OldMatrix, queue, "scale", scaleFactor);
            Future<HashMap<String, float[]>> submit = completionService.submit(worker);
            taskList.add(submit);
        }
        
        // collect results
        for (int t = 0; t < Threads; t++) {
            try {
                HashMap<String, float[]> hm = completionService.take().get();
                int r = 0;
                for (Map.Entry<String, float[]> entry : hm.entrySet()) {
                    String s = entry.getKey();
                    int i = Integer.valueOf(s);
                    int size = D - i;
                    float[] d = new float[size];
                    d = entry.getValue();
                    for (int j = 0; j < d.length; j++) {
                        int coord = j + i;
                        NewMatrix.setValueByEntry(d[j], i, coord);
                        r++;
                    }

                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        pool.shutdownNow();
        return NewMatrix;
    }
    
    public static GCNMatrix applySigmoid(GCNMatrix Similarity, float mu, float alpha, int Threads) {
        // prepwork
        int D = Similarity.getNumRows();
        String adjacency = "sigmoid";
        GCNMatrix Adjacency = new GCNMatrix(D, D);
        Adjacency = Operations.copyNames(Similarity.getRowNames(), Adjacency);
        Adjacency.setRowNames(Similarity.getRowNames());
        Adjacency.setAlternateRowNames(Similarity.getAlternateRowNames());
        Adjacency.setColumnNames(Similarity.getRowNames());
        // thread prep
        ExecutorService pool = Executors.newFixedThreadPool(Threads);
        ExecutorCompletionService<HashMap<String, float[]>> completionService = new ExecutorCompletionService<>(pool);
        List<Future<HashMap<String, float[]>>> taskList = new ArrayList<Future<HashMap<String, float[]>>>();
        ConcurrentLinkedQueue<String> queue = new ConcurrentLinkedQueue<String>();

        //queue prep
        for (int i = 0; i < D; i++) {
            String S = String.valueOf(i);
            queue.add(S);
        }
        int Number = D * D;
        System.err.println("Added "+ Number +" Tasks To Multithreaded Processing Engine");
        //add the tasks
        for (int i = 0; i < Threads; i++) {
            Callable<HashMap<String, float[]>> worker = new ConcurrentProcessing(Similarity, queue, adjacency, mu, alpha);
            Future<HashMap<String, float[]>> submit = completionService.submit(worker);
            taskList.add(submit);
        }
        
        // collect results
        for (int t = 0; t < Threads; t++) {
            try {
                HashMap<String, float[]> hm = completionService.take().get();
                int r = 0;
                for (Map.Entry<String, float[]> entry : hm.entrySet()) {
                    String s = entry.getKey();
                    int i = Integer.valueOf(s);
                    int size = D - i;
                    float[] d = new float[size];
                    d = entry.getValue();
                    for (int j = 0; j < d.length; j++) {
                        int coord = j + i;
                        Adjacency.setValueByEntry(d[j], i, coord);
                        r++;
                    }

                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        pool.shutdownNow();
        return Adjacency;
    }
    
    public static GCNMatrix calculateAdjacency(ExpressionFrame Expression, String corr, String adj, float mu, float alpha, int Threads,boolean signed) {

        // prepwork
        int D = Expression.getNumRows();
        GCNMatrix Adjacency = new GCNMatrix(D, D);
        Adjacency = Operations.copyNames(Expression.getRowNames(), Adjacency);
        switch (corr) {
            case "gini":
                Expression.calculateGiniSums(); // half of every gini coeff is a pre-calculable value
                break;
            case "pcc":
                Expression.calculateMeans(); // most of a pcc is pre-calculateable;
                break;
            case "spearman":
                break;
            default:
                explode();
                break;
        }

        // thread prep
        ExecutorService pool = Executors.newFixedThreadPool(Threads);
        ExecutorCompletionService<HashMap<String, float[]>> completionService = new ExecutorCompletionService<>(pool);
        List<Future<HashMap<String, float[]>>> taskList = new ArrayList<Future<HashMap<String, float[]>>>();
        ConcurrentLinkedQueue<String> queue = new ConcurrentLinkedQueue<String>();
		//System.err.println("Processing adjacency using " + Threads + " threads.");

        //queue prep
        for (int i = 0; i < D; i++) {
            String S = String.valueOf(i);
            queue.add(S);
        }
        int Number = D * D;
        System.err.println("Added "+ Number +" Tasks To Multithreaded Processing Engine");
        //add the tasks
        for (int i = 0; i < Threads; i++) {
            Callable<HashMap<String, float[]>> worker = new ConcurrentProcessing(Expression, queue, corr, "sigmoid", mu, alpha,signed);
            Future<HashMap<String, float[]>> submit = completionService.submit(worker);
            taskList.add(submit);
        }
        
        // collect results
        for (int t = 0; t < Threads; t++) {
            try {
                HashMap<String, float[]> hm = completionService.take().get();
                //System.err.println("obtained result for thread " + t);
                int r = 0;
                for (Map.Entry<String, float[]> entry : hm.entrySet()) {
                    String s = entry.getKey();
//					String[] S = s.split("-");
                    int i = Integer.valueOf(s);
                    int size = D - i;
                    float[] d = new float[size];
                    d = entry.getValue();
                    for (int j = 0; j < d.length; j++) {
                        int coord = j + i;
                        Adjacency.setValueByEntry(d[j], i, coord);
                        r++;
                    }

                }
                //System.err.println("Processed "+ r + " records");
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
            //System.err.println("Thread " + t + " complete.");
        }
        //System.err.println("Done.");
        pool.shutdownNow();
        return Adjacency;
    }

    public static GCNMatrix calculateDifferenceThreaded(GCNMatrix mat1, GCNMatrix mat2, int Threads) {

        // prepwork
        int D = mat1.getNumRows();
        GCNMatrix Difference = new GCNMatrix(D, D);
        Difference.setRowNames(mat1.getRowNames());
        Difference.setAlternateRowNames(mat2.getRowNames());
        Difference.setColumnNames(mat1.getRowNames());
        
        // thread prep
        ExecutorService pool = Executors.newFixedThreadPool(Threads);
        ExecutorCompletionService<HashMap<String, float[]>> completionService = new ExecutorCompletionService<>(pool);
        List<Future<HashMap<String, float[]>>> taskList = new ArrayList<Future<HashMap<String, float[]>>>();
        ConcurrentLinkedQueue<String> queue = new ConcurrentLinkedQueue<String>();
		//System.err.println("Processing adjacency using " + Threads + " threads.");

        //queue prep
        for (int i = 0; i < D; i++) {
            String S = String.valueOf(i);
            queue.add(S);
        }
        int Number = D * D;
        System.err.println("Added "+ Number +" Tasks To Multithreaded Processing Engine");
        //add the tasks
        for (int i = 0; i < Threads; i++) {
            Callable<HashMap<String, float[]>> worker = new ConcurrentProcessing(mat1,mat2, queue, "difference");
            Future<HashMap<String, float[]>> submit = completionService.submit(worker);
            taskList.add(submit);
        }
        
        // collect results
        for (int t = 0; t < Threads; t++) {
            try {
                HashMap<String, float[]> hm = completionService.take().get();
                //System.err.println("obtained result for thread " + t);
                int r = 0;
                for (Map.Entry<String, float[]> entry : hm.entrySet()) {
                    String s = entry.getKey();
                    int i = Integer.valueOf(s);
                    int size = D - i;
                    float[] d = new float[size];
                    d = entry.getValue();
                    for (int j = 0; j < d.length; j++) {
                        int coord = j + i;
                        Difference.setValueByEntry(d[j], i, coord);
                        //System.out.println(i+"\t"+j+"\t"+d);
                        r++;
                    }

                }
                //System.err.println("Processed "+ r + " records");
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
            //System.err.println("Thread " + t + " complete.");
        }
        //System.err.println("Done.");
        pool.shutdownNow();
        return Difference;
    }
    
    public static GCNMatrix calculateDifference(GCNMatrix mat1, GCNMatrix mat2) {
        int D = mat1.getNumRows();
        GCNMatrix Difference = new GCNMatrix(D, D);
        Difference.setRowNames(mat1.getRowNames());
        Difference.setAlternateRowNames(mat2.getRowNames());
        Difference.setColumnNames(mat1.getRowNames());
        for (int i = 0; i < D; i++) {
            for (int j = i; j < D; j++) {
                float v1 = mat1.getValueByEntry(i, j);
                float v2 = mat2.getValueByEntry(i, j);
                //if((v1 != 0) & (v2 != 0)){
                float d1 = v2 - v1;
                Difference.setValueByEntry(d1, i, j);
            }
        }
        return Difference;
    }

    public static void generateHistogramHM(GCNMatrix DataFrame, String pathOut, String Title, String Xlab, String Ylab, boolean print) {
        int H = DataFrame.getNumRows();
        int W = DataFrame.getNumColumns();
        DecimalFormat df = new DecimalFormat("#.##");
        df.setRoundingMode(RoundingMode.HALF_UP);
        TreeMap<Float, Integer> HMHistogram = new TreeMap<Float, Integer>();
        for (int i = 0; i < H; i++) {
            for (int j = 0; j < W; j++) {
                if (DataFrame.getValueByEntry(i, j) != 0) {
                    try {
                        Float V = (Float.valueOf(df.format(DataFrame.getValueByEntry(i, j))));
                        if (HMHistogram.containsKey(V)) {
                            Integer I = HMHistogram.get(V);
                            HMHistogram.put(V, I + 1);
                        } else {
                            HMHistogram.put(V, 1);
                        }
                    } catch (NumberFormatException ex) {
                        System.out.println("Obtain " + DataFrame.getValueByEntry(i, j) + " from matrix.");
                        System.exit(1);
                    }
                }
            }
        }
        XYSeriesCollection dataset = new XYSeriesCollection();
        XYSeries series = new XYSeries("Values");
        for (Map.Entry<Float, Integer> entry : HMHistogram.entrySet()) {
            Float A;
            A = entry.getKey();
            Integer value;
            value = entry.getValue();
            if (print == true) {
                System.out.println(A + "," + value);
            }
            series.add((double) A, (double) value, true);
        }
        dataset.addSeries(series);
        JFreeChart chart = ChartFactory.createXYLineChart(Title, Xlab, Ylab, dataset, PlotOrientation.VERTICAL,
                false, true, false);
        chart.setBackgroundPaint(Color.white);
        chart.setAntiAlias(true);
        final Plot plot = chart.getPlot();
        plot.setBackgroundPaint(Color.white);
        plot.setOutlinePaint(Color.black);
        try {
            ChartUtilities.saveChartAsJPEG(new File(pathOut), chart, 500, 300);
        } catch (IOException e) {
            System.err.println("Problem occurred creating chart.");
        }
    }


    public static float permuteData(ExpressionFrame expF1, ExpressionFrame expF2, int P, String out, String corr, float mu1, float mu2, float alpha1, float alpha2, int threads,boolean signed) {
        int s1 = expF1.getNumColumns();
        int s2 = expF2.getNumColumns();
        int R = expF1.getNumRows();
        int M = Math.min(expF1.getNumColumns(), expF2.getNumColumns());
        float pmu = (mu1+mu2)/2;
        float pa  = (alpha1+alpha2)/2;
        int S = s1 + s2;
        int s = M * 2;
        Integer[][] Sets = new Integer[P][];
        Sets = _getPermutations(s1, s2, P);
        float CUTOFF = 1.0f;
        ArrayList<TreeMap<Float, Integer>> Perms = new ArrayList<>();

        for (int p = 0; p < P; p++) {
            System.out.println("Permutation: " + p + " of " + P);
            ExpressionFrame pF1 = new ExpressionFrame(R, M);
            ExpressionFrame pF2 = new ExpressionFrame(R, M);
            for (int r = 0; r < R; r++) {
                float[] rF1 = expF1.getRowByIndex(r);
                float[] rF2 = expF2.getRowByIndex(r);
                float[] nR1 = new float[M];
                float[] nR2 = new float[M];
                for (int m = 0; m < M; m++) {
                    int ind = Sets[p][m];
                    if (ind < s1) {
                        nR1[m] = rF1[ind];
                    } else if (ind >= s1) {
                        ind = ind - s1;
                        nR1[m] = rF2[ind];
                    }
                }
                pF1.addRow(nR1);
                for (int m = M; m < s; m++) {
                    int Mind = m - M;
                    int ind = Sets[p][m];
                    if (ind < s1) {
                        nR2[Mind] = rF1[ind];
                    } else if (ind >= s1) {
                        ind = ind - s1;
                        nR2[Mind] = rF2[ind];
                    }
                }
                pF2.addRow(nR2);
            }
            System.out.println("Calculating Adjacencies on iteration " + p +"...");
            GCNMatrix CurrentMatrix1 = Operations.calculateAdjacency(pF1, corr, "sigmoid", pmu, pa, threads,signed);
            GCNMatrix CurrentMatrix2 = Operations.calculateAdjacency(pF2, corr, "sigmoid", pmu, pa, threads,signed);
            CurrentMatrix1.maskMatrix(0.01f);
            CurrentMatrix2.maskMatrix(0.01f);
            CurrentMatrix1.calculateKs();
            CurrentMatrix2.calculateKs();
            GCNMatrix Difference = Operations.calculateDifferenceThreaded(CurrentMatrix1, CurrentMatrix2,threads);
            TreeMap<Float, Integer> Distribution = Difference.generateDistribution();
            Perms.add(Distribution);
            System.out.println("Iteration "+ p +" complete.\n");
        }
        System.out.println("Calculating Observed Networks");
        System.out.println("Calculating Adjacencies...");
        GCNMatrix NetworkA = Operations.calculateAdjacency(expF1, corr, "sigmoid", mu1, alpha1, threads,signed);
        GCNMatrix NetworkB = Operations.calculateAdjacency(expF2, corr, "sigmoid", mu2, alpha2, threads,signed);
        NetworkA.maskMatrix(0.01f);
        NetworkB.maskMatrix(0.01f);
        NetworkA.calculateKs();
        NetworkB.calculateKs();
        GCNMatrix rDiff = Operations.calculateDifferenceThreaded(NetworkA, NetworkB,threads);
        TreeMap<Float, Integer> Real = rDiff.generateDistribution();
        String permutePathOut = out + "/PermutationDetails.tab";
        PrintWriter writer;
        DecimalFormat df = new DecimalFormat("#.######");
        df.setRoundingMode(RoundingMode.HALF_UP);
        try {
            writer = new PrintWriter(permutePathOut, "UTF-8");
            writer.println("Cutoff\tAverage False\tTrue\tFDR");
            for (float c = 0.0f; c < 2.0f; c += 0.01f) {
                float C = (Float.valueOf(df.format(c)));
                Double Total = 0.0d;
                for (int a = 0; a < Perms.size(); a++) {
                    for (Map.Entry<Float, Integer> entry : Perms.get(a).entrySet()) {
                        Float A;
                        A = entry.getKey();
                        Double value;
                        value = Double.valueOf(entry.getValue());
                        if (Math.abs(A) >= C) {
                            Total += value;
                        }
                    }
                }
                // Total holds all instances of adj value > C across all perms
                Double Average = Total / Perms.size();
                Double RealHits = 0.0d;
                for (Map.Entry<Float, Integer> entry : Real.entrySet()) {
                    Float A;
                    A = entry.getKey();
                    Double value;
                    value = Double.valueOf(entry.getValue());
                    if (Math.abs(A) >= C) {
                        RealHits += value;
                    }
                }
                double FDR = Average / RealHits;
                if(Double.isNaN(Average)) Average=0.0d;
                if(Double.isNaN(RealHits)) RealHits=0.0d;
                if(Double.isNaN(FDR)) FDR=0.0d;
                if(RealHits == 0.0d) FDR=1.0d;
                Average = (Double.valueOf(df.format(Average)));
                RealHits = (Double.valueOf(df.format(RealHits)));
                FDR = (Double.valueOf(df.format(FDR)));
                writer.println(C + "\t" + Average + "\t" + RealHits + "\t" + FDR);
                if (FDR <= 0.05f) {
                    if (CUTOFF == 1) {
                        CUTOFF = C;
                    } else {

                    }
                }

            }
            writer.close();
        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (UnsupportedEncodingException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return CUTOFF;
    }
    
    public static float[] permuteDataSigned(ExpressionFrame expF1, ExpressionFrame expF2, int P, String out, String corr, float neg_mu, float pos_mu, float neg_alpha, float pos_alpha, int threads,boolean signed) {
        int s1 = expF1.getNumColumns();
        int s2 = expF2.getNumColumns();
        int R = expF1.getNumRows();
        int M = Math.min(expF1.getNumColumns(), expF2.getNumColumns());
        int S = s1 + s2;
        int s = M * 2;
        Integer[][] Sets = new Integer[P][];
        Sets = _getPermutations(s1, s2, P);
        float CUTOFF[] = new float[2];
        CUTOFF[0]=1.0f;
        CUTOFF[1]=1.0f;
        ArrayList<TreeMap<Float, Integer>> NegativePerms = new ArrayList<>();
        ArrayList<TreeMap<Float, Integer>> PositivePerms = new ArrayList<>();

        for (int p = 0; p < P; p++) {
            System.out.println("Permutation: " + p + " of " + P);
            ExpressionFrame pF1 = new ExpressionFrame(R, M);
            ExpressionFrame pF2 = new ExpressionFrame(R, M);
            for (int r = 0; r < R; r++) {
                float[] rF1 = expF1.getRowByIndex(r);
                float[] rF2 = expF2.getRowByIndex(r);
                float[] nR1 = new float[M];
                float[] nR2 = new float[M];
                for (int m = 0; m < M; m++) {
                    int ind = Sets[p][m];
                    if (ind < s1) {
                        nR1[m] = rF1[ind];
                    } else if (ind >= s1) {
                        ind = ind - s1;
                        nR1[m] = rF2[ind];
                    }
                }
                pF1.addRow(nR1);
                for (int m = M; m < s; m++) {
                    int Mind = m - M;
                    int ind = Sets[p][m];
                    if (ind < s1) {
                        nR2[Mind] = rF1[ind];
                    } else if (ind >= s1) {
                        ind = ind - s1;
                        nR2[Mind] = rF2[ind];
                    }
                }
                pF2.addRow(nR2);
            }
            System.out.println("Calculating Adjacencies on iteration " + p +"...");
            GCNMatrix CurrentMatrix1 = Operations.calculateAdjacency(pF1, corr, "sigmoid", 0.0f, 0.0f, threads,signed);
            GCNMatrix CurrentMatrix2 = Operations.calculateAdjacency(pF2, corr, "sigmoid", 0.0f, 0.0f, threads,signed);
            CurrentMatrix1.maskMatrix(0.01f);
            CurrentMatrix2.maskMatrix(0.01f);
            CurrentMatrix1.calculateKs();
            CurrentMatrix2.calculateKs();
            
            GCNMatrix Difference = Operations.calculateDifferenceThreaded(CurrentMatrix1, CurrentMatrix2,threads);
            Difference = Operations.findElasticity(Difference,0.0f,CurrentMatrix1, CurrentMatrix2,"negative");
            TreeMap<Float, Integer> Distribution = Difference.generateSignedDistribution();
            NegativePerms.add(Distribution);
            
            Difference = null;
            Difference = Operations.calculateDifferenceThreaded(CurrentMatrix1, CurrentMatrix2,threads);
            Difference = Operations.findElasticity(Difference,0.0f,CurrentMatrix1, CurrentMatrix2,"positive");
            Distribution = Difference.generateSignedDistribution();
            PositivePerms.add(Distribution);
            
            int disp_P = p + 1;
            System.out.println("Iteration "+ disp_P +" complete.\n");
        }
        System.out.println("Calculating Observed Networks");
        System.out.println("Calculating Adjacencies...");
        GCNMatrix NetworkA = Operations.calculateAdjacency(expF1, corr, "sigmoid", 0.0f, 0.0f, threads,signed);
        GCNMatrix NetworkB = Operations.calculateAdjacency(expF2, corr, "sigmoid", 0.0f, 0.0f, threads,signed);
        NetworkA.maskMatrix(0.01f);
        NetworkB.maskMatrix(0.01f);
        NetworkA.calculateKs();
        NetworkB.calculateKs();
        GCNMatrix rDiff = Operations.calculateDifferenceThreaded(NetworkA, NetworkB, threads);
        rDiff = Operations.findElasticity(rDiff,0.0f,NetworkA, NetworkB,"negative");
        TreeMap<Float, Integer> RealNeg = rDiff.generateSignedDistribution();
        
        rDiff=null;
        rDiff = Operations.calculateDifferenceThreaded(NetworkA, NetworkB, threads);
        rDiff = Operations.findElasticity(rDiff,0.0f,NetworkA, NetworkB,"positive");
        TreeMap<Float, Integer> RealPos = rDiff.generateSignedDistribution();
        
        String permutePathOut_Adj = out + "/FDR_Calculations_Adjacency.tab";
        String permutePathOut_Sim = out + "/FDR_Calculations_Similarity.tab";
        PrintWriter writer_Adj;
        PrintWriter writer_Sim;
        DecimalFormat df = new DecimalFormat("#.######");
        df.setRoundingMode(RoundingMode.HALF_UP);
        try {
            writer_Adj = new PrintWriter(permutePathOut_Adj, "UTF-8");
            writer_Adj.println("Cutoff\tAverage False\tTrue\tFDR");
            writer_Sim = new PrintWriter(permutePathOut_Sim, "UTF-8");
            writer_Sim.println("Cutoff\tAverage False\tTrue\tFDR");
            for (float c = -2.0f; c < 0.0f; c += 0.01f) {
                float C = (Float.valueOf(df.format(c)));
                Double Total = 0.0d;
                for (int a = 0; a < NegativePerms.size(); a++) {
                    for (Map.Entry<Float, Integer> entry : NegativePerms.get(a).entrySet()) {
                        Float A = entry.getKey();
                        Double value = Double.valueOf(entry.getValue());
                        if (A <= C){
                            Total += value;
                        }
                    }
                }
                // Total holds all instances of adj value > C across all perms
                Double Average = Total / NegativePerms.size();
                Double RealHits = 0.0d;
                for (Map.Entry<Float, Integer> entry : RealNeg.entrySet()) {
                    Float A = entry.getKey();
                    Double value = Double.valueOf(entry.getValue());
                    if (A <= C) {
                        RealHits += value;
                    }
                }
                double FDR = Average / RealHits;
                
                if(Double.isNaN(Average)) Average=0.0d;
                if(Double.isNaN(RealHits)) RealHits=0.0d;
                if(Double.isNaN(FDR)) FDR=0.0d;
                if(RealHits == 0.0d) FDR=1.0d;
                Average = (Double.valueOf(df.format(Average)));
                RealHits = (Double.valueOf(df.format(RealHits)));
                FDR = (Double.valueOf(df.format(FDR)));
                writer_Sim.println(C + "\t" + Average + "\t" + RealHits + "\t" + FDR);
                float A=C/2.0f;
                A = _getSigmoid(A,neg_mu,neg_alpha);
                A = A * 2.0f;
                if(Math.abs(A) >= 0.01f){
                    writer_Adj.println(A + "\t" + Average + "\t" + RealHits + "\t" + FDR);
                }
                if ((FDR >= 0.05f) && (FDR != 1.0f)) {
                    if (CUTOFF[0] == 1.0f) {
                        CUTOFF[0] = A-0.01f;
                    } else {

                    }
                }
            }
            for (float c = 0.0f; c < 2.0f; c += 0.01f) {
                float C = (Float.valueOf(df.format(c)));
                Double Total = 0.0d;
                for (int a = 0; a < PositivePerms.size(); a++) {
                    for (Map.Entry<Float, Integer> entry : PositivePerms.get(a).entrySet()) {
                        Float A = entry.getKey();
                        Double value = Double.valueOf(entry.getValue());
                        if (A >= C){
                            Total += value;
                        }
                    }
                }
                // Total holds all instances of adj value > C across all perms
                Double Average = Total / PositivePerms.size();
                Double RealHits = 0.0d;
                for (Map.Entry<Float, Integer> entry : RealPos.entrySet()) {
                    Float A = entry.getKey();
                    Double value = Double.valueOf(entry.getValue());
                    if (A >= C) {
                        RealHits += value;
                    }
                }
                double FDR = Average / RealHits;
                
                if(Double.isNaN(Average)) Average=0.0d;
                if(Double.isNaN(RealHits)) RealHits=0.0d;
                if(Double.isNaN(FDR)) FDR=0.0d;
                if(RealHits == 0.0d) FDR=1.0d;
                Average = (Double.valueOf(df.format(Average)));
                RealHits = (Double.valueOf(df.format(RealHits)));
                FDR = (Double.valueOf(df.format(FDR)));
                writer_Sim.println(C + "\t" + Average + "\t" + RealHits + "\t" + FDR);
                float A=C/2.0f;
                A = _getSigmoid(A,pos_mu,pos_alpha);
                A = A * 2.0f;
                if(Math.abs(A) >= 0.01f){
                    writer_Adj.println(A + "\t" + Average + "\t" + RealHits + "\t" + FDR);
                }
                if (FDR <= 0.05f) {
                    if (CUTOFF[1] == 1.0f) {
                        CUTOFF[1] = A;
                    } else {

                    }
                }
            }
            
            writer_Sim.close();
            writer_Adj.close();
        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (UnsupportedEncodingException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return CUTOFF;
    }
    private static float _getSigmoid(float V, float M, float A) {
        if((A == 0.0f) && (M == 0.0f)){
            return V; // passthrough - must be notated in manual as A=0,M=0 =/= sigmoid(V)=V
        }else{
            float s = (float) (1.0f / (1 + Math.exp(A * -1 * (Math.abs(V) - M))));
            if(V<0.0f){
                s=s*-1;
            }
            return s;
        }
        
    }

    public static float[] determineCutoffSF(ExpressionFrame expF1, ExpressionFrame expF2, String out, String corr, float mu1, float mu2, float alpha1, float alpha2, int threads) throws FileNotFoundException, UnsupportedEncodingException {
        boolean signed=true;
        GCNMatrix CurrentMatrix1 = Operations.calculateAdjacency(expF1, corr, "sigmoid", mu1, alpha1, threads,signed);
        GCNMatrix CurrentMatrix2 = Operations.calculateAdjacency(expF2, corr, "sigmoid", mu2, alpha2, threads,signed);
        CurrentMatrix1.calculateKs();
        CurrentMatrix2.calculateKs();
        
        float[] result = new float[2];
        
        PrintWriter writer;
        try{
            DecimalFormat df = new DecimalFormat("#.###");
            df.setRoundingMode(RoundingMode.HALF_UP);
            String permutePathOut = out + "/NegativePermutationDetails.tab";
            writer= new PrintWriter(permutePathOut, "UTF-8");
            writer.println("Cutoff\tR-squared\tSlope\tMean Connectivity");    
            System.err.println("Determining Cutoff for Negative Elasticity.");
            GCNMatrix Difference = Operations.calculateDifferenceThreaded(CurrentMatrix1, CurrentMatrix2,threads);
            for(float cutoff=0;cutoff<2.0f;cutoff+=0.01){
                cutoff = (Float.valueOf(df.format(cutoff)));
                float this_cutoff = cutoff *-1.0f;
                //Difference.maskAbove(this_cutoff);
                Difference = Operations.findElasticity(Difference,this_cutoff,CurrentMatrix1,CurrentMatrix2,"negative");
                Difference.calculateKs();
                double[] Return = Difference.determineScaleFreeCritereon();
                double RSquared = Return[0];
                double Slope = Return[1];
                float mean = Difference.getMeanK();
                if(Double.isNaN(RSquared)) RSquared=0.0d;
                RSquared=(Double.valueOf(df.format(RSquared)));
                mean = (Float.valueOf(df.format(mean)));
                writer.println(this_cutoff + "\t" + RSquared + "\t" + Slope + "\t" + mean);
                if(result[0] == 0){
                    if((RSquared>=0.8)&&(Slope < -0.8)) result[0]=cutoff;
                }
            }
            writer.close();
            permutePathOut = out + "/PositivePermutationDetails.tab";
            System.err.println("Determining Cutoff for Positive Elasticity.");
            Difference = Operations.calculateDifferenceThreaded(CurrentMatrix1, CurrentMatrix2,threads);
            writer= new PrintWriter(permutePathOut, "UTF-8");
            writer.println("Cutoff\tR-squared\tSlope\tMean Connectivity");    
            for(float cutoff=0;cutoff<2.0f;cutoff+=0.01){
                cutoff = (Float.valueOf(df.format(cutoff)));
                //Difference.maskBelow(cutoff);
                Difference = Operations.findElasticity(Difference,cutoff,CurrentMatrix1,CurrentMatrix2,"positive");
                Difference.calculateKs();
                double[] Return = Difference.determineScaleFreeCritereon();
                double RSquared = Return[0];
                double Slope = Return[1];
                float mean = Difference.getMeanK();
                if(Double.isNaN(RSquared)) RSquared=0.0d;
                RSquared=(Double.valueOf(df.format(RSquared)));
                mean = (Float.valueOf(df.format(mean)));
                writer.println(cutoff + "\t" + RSquared + "\t" + Slope + "\t" + mean);
                if(result[1] == 0){
                    if((RSquared>=0.8)&&(Slope < -0.8)) result[1]=cutoff;
                }
            }
            writer.close();
        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (UnsupportedEncodingException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return result;
    }

    
    public static GCNMatrix findElasticity (GCNMatrix Difference, float Cut, GCNMatrix Net1, GCNMatrix Net2, String type){
        int this_N = Difference.getN();
        GCNMatrix Elasticity = new GCNMatrix(this_N, this_N);
        Elasticity.setRowNames(Net1.getRowNames());
        Elasticity.setAlternateRowNames(Net2.getRowNames());
        Elasticity.setColumnNames(Net1.getRowNames());
        for(int i=0;i<this_N;i++){
            for(int j=i;j<this_N;j++){
                float d = Math.abs(Difference.getValueByEntry(i,j));
                if(d < Math.abs(Cut)) continue;
                float v1 = Net1.getValueByEntry(i,j);
                float v2 = Net2.getValueByEntry(i,j);
                if("positive".equals(type)){
                    if(Math.abs(v2)<Math.abs(v1)) continue;
                    Elasticity.setValueByEntry(d, i, j);
                }
                if("negative".equals(type)){
                    if(Math.abs(v2)>Math.abs(v1)) continue;
                    d = d*-1.0f;
                    Elasticity.setValueByEntry(d, i, j);
                }
                
            }
        }
        return Elasticity;
    }
    
    private static Integer[][] _getPermutations(int s1, int s2, int p) {
        Integer[][] Sets = new Integer[p][];
        int m = Math.min(s1, s2);
        int S = s1 + s2;
        int s = m * 2;
        Integer ind[] = new Integer[S];
        for (int i = 0; i < S; i++) {
            ind[i] = i;
        }
        for (int i = 0; i < p; i++) {
            Integer shuffled[] = fisherYates(ind);
            Integer[] perms = Arrays.copyOfRange(shuffled, 0, s);
            Sets[i] = perms;
        }
        return Sets;
    }

    /// shuffle some primitives
    private static Integer[] fisherYates(Integer[] array) {
        Random rnd = new Random();
        for (int i = array.length - 1; i > 0; i--) {
            int index = rnd.nextInt(i + 1);
            // Simple swap
            int a = array[index];
            array[index] = array[i];
            array[i] = a;
        }
        return array;
    }

    public static void createDirectory(String dir) {
        File theDir = new File(dir);
        if (!theDir.exists()) {
            System.out.println("creating directory: " + dir);
            try {
                theDir.mkdir();
            } catch (SecurityException se) {
                //TODO handle it
            }
        }
    }

    private static void _tempPrintPermsToFile(float[] rTOMS, float[][] TOMpermutations, String outBase, String[] names) {
        try {
            String path = outBase + "/Selfwise.pergene.TOM.tab";
            String Sep = "\t";
            PrintWriter writer = new PrintWriter(path, "UTF-8");
            for (int g = 0; g < rTOMS.length; g++) {
                writer.print(names[g] + Sep + rTOMS[g]);
                for (int p = 0; p < TOMpermutations.length; p++) {
                    writer.print(Sep + TOMpermutations[p][g]);
                }
                writer.print("\n");
            }

            writer.close();
        } catch (Exception e) {
            // 
        }
    }

}
