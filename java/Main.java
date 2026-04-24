import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Main {
    public static void main(String[] args) throws IOException {
        int N = 200;
        if (args.length > 0) N = Integer.parseInt(args[0]);
        boolean timing = args.length > 1 && args[1].equals("timing");

        double RADIUS = 1;
        double MASS = 1.0;
        double INITIAL_SPEED = 1;
        double MAX_TIME = timing ? 5.0 : 1000.0;

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
                    continue;
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
        if (timing) {
            system.simulateTiming(MAX_TIME);
        } else {
            system.simulate(MAX_TIME, 1);
        }
    }
}
