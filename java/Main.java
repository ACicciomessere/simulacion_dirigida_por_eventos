import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Main {
    public static void main(String[] args) throws IOException {
        // --- Simulation Parameters ---
        if (args.length == 0) {
            System.err.println("Usage: java Main <N> [timing]");
            System.exit(1);
        }

        int N = Integer.parseInt(args[0]); // Número de partículas
        boolean timingMode = args.length > 1 && args[1].equals("timing");

        double RADIUS = 1; // Radio de la partícula
        double MASS = 1.0; // Masa de la partícula
        double INITIAL_SPEED = 1; // Magnitud de la velocidad inicial
        double MAX_TIME = timingMode ? 500.0 : 1000.0; // Timing mode: 5s, else full simulation

        double R_OUTER = 40;
        double R_INNER = 1;
        double CX = 0.0;
        double CY = 0.0;

        // INFO: 1. se definen las posiciones y velocidades iniciales, los radios y
        // tamaño de la caja
        List<Particle> particles = new ArrayList<>();
        Random rand = new Random();

        for (int i = 0; i < N; i++) {
            double x, y;
            boolean overlap;

            do {
                overlap = false;

                x = CX - R_OUTER + rand.nextDouble() * (2 * R_OUTER);
                y = CY - R_OUTER + rand.nextDouble() * (2 * R_OUTER);

                // 1.1 que cada particula nueva x no se superponga con ninguna de las existentes
                // ni con las paredes

                // paredes
                double distToCenter = Math.sqrt((x - CX) * (x - CX) + (y - CY) * (y - CY));

                if (distToCenter > (R_OUTER - RADIUS) ||
                        distToCenter < (R_INNER + RADIUS)) {
                    overlap = true;
                }

                // entre partículas
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

        // --- Ejecutar Simulación ---
        Collisions system = new Collisions(particles, R_OUTER, R_INNER, CX, CY);
        system.simulate(MAX_TIME, timingMode ? 1000 : 1);
    }
}
