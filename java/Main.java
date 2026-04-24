import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Main {
    public static void main(String[] args) throws IOException {
        // --- Simulation Parameters ---
        int N = args.length > 0 ? Integer.parseInt(args[0]) : 200;
        // Output file: pass as second arg, or default to ../outputs/sim_circular/output.txt
        String outputPath = args.length > 1 ? args[1] : "../outputs/sim_circular/output.txt";

        double RADIUS = 1;
        double MASS = 1.0;
        double INITIAL_SPEED = 1;
        double MAX_TIME = 1000.0;
        int N_EVENTS_PER_PRINT = 1; // snapshot at every bounce event (required by analyze.py)

        double R_OUTER = 40;
        double R_INNER = 1;
        double CX = 0.0;
        double CY = 0.0;

        List<Particle> particles = new ArrayList<>();
        Random rand = new Random();

        for (int i = 0; i < N; i++) {
            double x, y;
            boolean overlap;

            do {
                overlap = false;
                x = CX - R_OUTER + rand.nextDouble() * (2 * R_OUTER);
                y = CY - R_OUTER + rand.nextDouble() * (2 * R_OUTER);

                double distToCenter = Math.sqrt((x - CX) * (x - CX) + (y - CY) * (y - CY));
                if (distToCenter > (R_OUTER - RADIUS) || distToCenter < (R_INNER + RADIUS)) {
                    overlap = true;
                }

                for (Particle existing : particles) {
                    double distSq = Math.pow(x - existing.x, 2) + Math.pow(y - existing.y, 2);
                    if (distSq < Math.pow(2 * RADIUS, 2)) {
                        overlap = true;
                        break;
                    }
                }
            } while (overlap);

            double angle = 2 * Math.PI * rand.nextDouble();
            double vx = INITIAL_SPEED * Math.cos(angle);
            double vy = INITIAL_SPEED * Math.sin(angle);

            particles.add(new Particle(x, y, vx, vy, RADIUS, MASS));
        }

        Collisions system = new Collisions(particles, R_OUTER, R_INNER, CX, CY);
        system.simulate(MAX_TIME, N_EVENTS_PER_PRINT, outputPath);
    }
}
