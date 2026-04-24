import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Locale;
import java.util.PriorityQueue;

public class Collisions {
    private PriorityQueue<Event> pq;
    private double simTime;
    private final List<Particle> particles;

    private final double R_OUTER;
    private final double R_INNER;
    private final double CX;
    private final double CY;

    public Collisions(List<Particle> particles, double rOuter, double rInner, double cx, double cy) {
        this.particles = particles;
        this.R_OUTER = rOuter;
        this.R_INNER = rInner;
        this.CX = cx;
        this.CY = cy;
        this.simTime = 0.0;
        this.pq = new PriorityQueue<>();
    }

    private void predict(Particle p) {
        if (p == null) return;

        double tOuter = p.timeToHitCircle(CX, CY, R_OUTER, false);
        if (tOuter > 1e-12) {
            pq.add(new Event(simTime + tOuter, p, null));
        }

        double tInner = p.timeToHitCircle(CX, CY, R_INNER, true);
        if (tInner > 1e-12) {
            pq.add(new Event(simTime + tInner, p, p));
        }
    }

    // Full simulation with file output
    public void simulate(double maxTime, int nEventsPerPrint) throws IOException {
        for (Particle p : particles) predict(p);
        pq.add(new Event(0, null, null));

        String simPath = "../outputs/sim_circular";
        File folder = new File(simPath);
        if (!folder.exists()) folder.mkdirs();

        try (FileWriter outputWriter = new FileWriter(simPath + "/output.txt")) {
            outputWriter.write("# t\n");
            outputWriter.write("# x y vx vy state\n");

            int countEvents = 0;
            while (!pq.isEmpty()) {
                Event event = pq.poll();
                if (!event.isValid()) continue;
                if (event.time < simTime - 1e-12) continue;

                for (Particle p : particles) p.move(event.time - simTime);
                simTime = event.time;
                countEvents++;

                if (countEvents % nEventsPerPrint == 0) {
                    outputWriter.write(String.format(Locale.US, "%.6f\n", simTime));
                    for (Particle p : particles) {
                        outputWriter.write(String.format(Locale.US,
                                "%.6f %.6f %.6f %.6f %d\n",
                                p.x, p.y, p.vx, p.vy, p.state));
                    }
                }

                if (simTime > maxTime) break;

                resolveEvent(event);
            }
        }
        System.out.println("Simulation finished.");
    }

    // Timing mode: no file I/O, fixed tf=5s
    public void simulateTiming(double maxTime) {
        for (Particle p : particles) predict(p);
        pq.add(new Event(0, null, null));

        while (!pq.isEmpty()) {
            Event event = pq.poll();
            if (!event.isValid()) continue;
            if (event.time < simTime - 1e-12) continue;

            for (Particle p : particles) p.move(event.time - simTime);
            simTime = event.time;

            if (simTime > maxTime) break;

            resolveEvent(event);
        }
        System.out.println("Timing finished.");
    }

    private void resolveEvent(Event event) {
        Particle a = event.a;
        Particle b = event.b;

        if (a != null && b != null) {
            if (a == b) {
                // Colisión con obstáculo interno: fresca → usada
                a.bounceOffCircle(CX, CY);
                a.state = 1;
            } else {
                // Colisión partícula-partícula
                a.bounceOff(b);
            }
        } else if (a != null) {
            // Colisión con pared externa: usada → fresca
            a.bounceOffCircle(CX, CY);
            if (a.state == 1) a.state = 0;
        }

        predict(a);
        if (b != null && b != a) predict(b);
    }
}
