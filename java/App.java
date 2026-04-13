import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class App {
    public static void main(String[] args) {

        if (args.length < 7) {
            System.err.println(
                    "HELP (File Mode): java -cp src App <StaticPath> <DynamicPath> <M> <rc> <periodic (true/false)> <iterations> <r0>");
            System.err.println(
                    "HELP (Random Mode): java -cp src App <N> <M> <rc> <periodic (true/false)> <iterations> <r0>");
            System.exit(1);
        }

        // recibimos los parametros desde la terminal, al correr run.sh
        int M = Integer.parseInt(args[1]);
        double rc = Double.parseDouble(args[2]);
        boolean periodic = Boolean.parseBoolean(args[3]);
        int iterations = Integer.parseInt(args[4]);
        
        double r0 = Double.parseDouble(args[args.length - 1]);
        double r_inner = 1.0; // Inner circle radius
        boolean circleLeader = false;

        if (args.length > 7) {
            circleLeader = Boolean.parseBoolean(args[5]);
        }

        ArrayList<Particle> particles = new ArrayList<>();
        int N = 0;
        double L = 2.0 * r0;  // For compatibility, but circles are centered at origin

        boolean isRandomMode = false;
        try {
            // tmb parametros desde la terminal, en este caso serian números fijos
            N = Integer.parseInt(args[0]);
            isRandomMode = true;
        } catch (NumberFormatException e) {
            isRandomMode = false;
        }

        if (isRandomMode) {
            System.out.println("Running in Random Mode...");
            // los limites marcados en el tp1
            double r_min = 0.23;
            double r_max = 0.26;
            double property = 1.0;
            double theta = 0.0;
            java.util.Random rand = new java.util.Random();

            // le inventamos posiciones y radios a las N partículas q pedimos, chequeando
            // que no se superpongan entre si.
            for (int i = 0; i < N; i++) {
                boolean overlaps;
                double rx, ry, radius;
                int attempts = 0;

                do {
                    overlaps = false;
                    radius = r_min + rand.nextDouble() * (r_max - r_min);
                    
                    double max_r = r0 - radius;
                    if (max_r < 0) max_r = 0;
                    double r_dist = Math.sqrt(rand.nextDouble()) * max_r;
                    double angle = rand.nextDouble() * 2 * Math.PI;
                    
                    rx = r_dist * Math.cos(angle);
                    ry = r_dist * Math.sin(angle);
                    
                    theta = rand.nextDouble() * 2 * Math.PI;

                    for (Particle p : particles) {
                        double dx = Math.abs(rx - p.getX());
                        double dy = Math.abs(ry - p.getY());

                        if (periodic) {
                            if (dx > L / 2)
                                dx = L - dx;
                            if (dy > L / 2)
                                dy = L - dy;
                        }

                        double minDistance = radius + p.getRadius();
                        if (Math.sqrt(dx * dx + dy * dy) < minDistance) {
                            overlaps = true;
                            break;
                        }
                    }
                    attempts++;
                } while (overlaps && attempts < 10000);

                particles.add(new Particle(i, rx, ry, theta, radius, property, false));
                
            }
        } else {
            System.out.println("Running in File Mode...");
            // tmb desde la terminal, en este caso serian archivos
            String staticPath = args[0];
            String dynamicPath = args[1];

            try {
                Scanner staticScanner = new Scanner(new File(staticPath));
                Scanner dynamicScanner = new Scanner(new File(dynamicPath));

                N = Integer.parseInt(staticScanner.nextLine().trim());

                for (int i = 0; i < N; i++) {
                    String[] staticParts = staticScanner.nextLine().trim().split("\\s+");
                    String[] dynamicParts = dynamicScanner.nextLine().trim().split("\\s+");

                    double radius = Double.parseDouble(staticParts[0]);
                    double property = Double.parseDouble(staticParts[1]);

                    double rx = Double.parseDouble(dynamicParts[0]);
                    double ry = Double.parseDouble(dynamicParts[1]);
                    double theta = Double.parseDouble(dynamicParts[2]);
                    // Normalize to center at origin for circular boundary
                    rx = rx - L / 2.0;
                    ry = ry - L / 2.0;

                    // en este caso no se hacen todos los chequeos como en el anterior porq se asume
                    // q los archivos son correctos TODO
                    particles.add(new Particle(i, rx, ry, theta, radius, property, false));
                }

                staticScanner.close();
                dynamicScanner.close();

            } catch (FileNotFoundException e) {
                System.err.println("Error reading input files: " + e.getMessage());
                System.exit(1);
            }
        }

        System.out.println("Successfully parsed " + N + " particles.");

        double cellSize = 2.0 * r0 / M;
        if (cellSize <= rc) {
            //TODO, esta condición se mantiene al anterior TP??
            System.err.println("Error: The condition L/M > rc is not met. (L/M=" + cellSize + ", rc=" + rc + ")");
            System.err.println("Execution stopped to prevent loss of neighbors.");
            System.exit(1);
        }

        // Ejecuta Cell Index Method para calcular los vecinos
        automataCelular(particles, L, M, rc, periodic, iterations, circleLeader, r0, r_inner);
    }

    public static void automataCelular(ArrayList<Particle> particles, double L, int M, double rc, boolean periodic,
            int iterations, boolean circleLeader, double r0, double r_inner) {
        
        java.util.Random rand = new java.util.Random();
        System.out.println("Starting simulation for " + iterations + " iterations...");
        double ETA = 0.1; //TODO eta no deberia existir

        exportFrame(0, particles, L); // Exportamos el estado inicial verdadero

        for (int t = 1; t <= iterations; t++) {
            cellIndexMethod(particles, L, M, rc, periodic, r0);

            for (Particle p : particles) {
                p.calculateNextTheta(ETA, rand);
            }

            for (Particle p : particles) {
                p.updateTheta();
                p.updatePosition(r0, r_inner, periodic);
            }

            exportFrame(t, particles, L);
        }
        System.out.println("Simulation finished.");
    }

    public static double getDistance(Particle p1, Particle p2, double L, boolean periodic) {
        double dx = Math.abs(p1.getX() - p2.getX());
        double dy = Math.abs(p1.getY() - p2.getY());

        if (periodic) {
            if (dx > L / 2)
                dx = L - dx;
            if (dy > L / 2)
                dy = L - dy;
        }

        return Math.sqrt(dx * dx + dy * dy) - p1.getRadius() - p2.getRadius();
    }

    // se usa en benchmark.java para comparar tiempos de ejecución entre ambos
    // métodos
    public static long bruteForce(ArrayList<Particle> particles, double L, double rc, boolean periodic) {
        for (Particle p : particles)
            p.clearNeighbours();

        long start = System.nanoTime();
        for (int i = 0; i < particles.size(); i++) {
            Particle p1 = particles.get(i);
            for (int j = i + 1; j < particles.size(); j++) {
                Particle p2 = particles.get(j);
                if (getDistance(p1, p2, L, periodic) <= rc) {
                    p1.addNeighbour(p2);
                    p2.addNeighbour(p1);
                }
            }
        }
        long end = System.nanoTime();
        long diff = end - start;
        System.out.println("Brute Force Execution Time: " + (diff / 1_000_000.0) + " ms");
        return diff;
    }

    public static long cellIndexMethod(ArrayList<Particle> particles, double L, int M, double rc, boolean periodic, double r0) {
        for (Particle p : particles)
            p.clearNeighbours();

        long start = System.nanoTime();
        ArrayList<ArrayList<Particle>> cells = new ArrayList<>(M * M);
        for (int i = 0; i < M * M; i++) {
            cells.add(new ArrayList<>());
        }

        double r_max = 0;
        for (Particle p : particles) {
            if (p.getRadius() > r_max)
                r_max = p.getRadius();
        }

        double cellSize = L / M;
        // Calculamos cuántas celdas hacia los lados debemos iterar
        int searchRange = (int) Math.ceil((rc + 2 * r_max) / cellSize);

        for (Particle p : particles) {
            // With this:
            int cellX = (int) Math.floor((p.getX() + r0) / cellSize);
            int cellY = (int) Math.floor((p.getY() + r0) / cellSize);

            if (cellX == M)
                cellX = M - 1;
            if (cellY == M)
                cellY = M - 1;

            cells.get(cellY * M + cellX).add(p);
        }

        for (int cy = 0; cy < M; cy++) {
            for (int cx = 0; cx < M; cx++) {
                List<Particle> cell1 = cells.get(cy * M + cx);
                if (cell1.isEmpty())
                    continue;

                for (int dy = -searchRange; dy <= searchRange; dy++) {
                    for (int dx = -searchRange; dx <= searchRange; dx++) {
                        int nx = cx + dx;
                        int ny = cy + dy;

                        if (!periodic) {
                            if (nx < 0 || nx >= M || ny < 0 || ny >= M)
                                continue;
                        } else {
                            nx = (nx % M + M) % M;
                            ny = (ny % M + M) % M;
                        }

                        List<Particle> cell2 = cells.get(ny * M + nx);
                        if (cell2.isEmpty())
                            continue;

                        for (Particle p1 : cell1) {
                            for (Particle p2 : cell2) {
                                // Evita comparar dos veces la misma pareja, y evita comparar una partícula
                                // consigo misma
                                if (p1.getId() < p2.getId()) {
                                    if (getDistance(p1, p2, L, periodic) <= rc) {
                                        // ahorramos iteraciones al agregar a ambos vecinos entre si al mismo tiempo
                                        p1.addNeighbour(p2);
                                        p2.addNeighbour(p1);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        long end = System.nanoTime();
        long diff = end - start;
        System.out.println("Cell Index Method Execution Time: " + (diff / 1_000_000.0) + " ms");
        return diff;
    }

    private static void exportFrame(int t, ArrayList<Particle> particles, double L) {
        try (FileWriter writer = new FileWriter("../particles_frames.txt", true)) {
            writer.write(particles.size() + "\n");
            writer.write("Frame " + t + "\n");
            for (Particle p : particles) {
                int leaderFlag = p.isLeader() ? 1 : 0;
                writer.write(p.getId() + " " + p.getX() + " " + p.getY() + " "
                        + Math.cos(p.getTheta()) * Particle.VELOCITY + " " + Math.sin(p.getTheta()) * Particle.VELOCITY
                        + " " + p.getRadius() + " " + leaderFlag + "\n");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}