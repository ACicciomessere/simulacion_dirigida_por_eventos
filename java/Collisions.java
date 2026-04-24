import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.PriorityQueue;

public class Collisions {
    private PriorityQueue<Event> pq;
    private double simTime;
    private final List<Particle> particles;

    // Geometría circular
    private final double R_OUTER;
    private final double R_INNER;
    private final double CX;
    private final double CY;

    // Parámetro para limpiar eventos inválidos
    private static final int CLEANUP_THRESHOLD = 10000; // Limpiar cuando la cola exceda esto
    private int eventsSinceCleanup = 0;

    public Collisions(List<Particle> particles, double rOuter, double rInner, double cx, double cy) {
        this.particles = particles;
        this.R_OUTER = rOuter;
        this.R_INNER = rInner;
        this.CX = cx;
        this.CY = cy;

        this.simTime = 0.0;
        this.pq = new PriorityQueue<>();
    }

    private void cleanupQueue() {
        // Reconstruir la cola quitando eventos inválidos
        List<Event> validEvents = new ArrayList<>();
        while (!pq.isEmpty()) {
            Event e = pq.poll();
            if (e.isValid()) {
                validEvents.add(e);
            }
        }
        pq = new PriorityQueue<>(validEvents);
        eventsSinceCleanup = 0;
    }

    private void predict(Particle p) {
        if (p == null)
            return;

        // 🔹 pared externa
        double tOuter = p.timeToHitCircle(CX, CY, R_OUTER, false);
        if (tOuter > 1e-12) {
            pq.add(new Event(simTime + tOuter, p, null));
        }

        // 🔹 círculo interno
        double tInner = p.timeToHitCircle(CX, CY, R_INNER, true);
        if (tInner > 1e-12) {
            pq.add(new Event(simTime + tInner, p, p)); // usamos (a==b) como flag
        }

        // 🔹 otras partículas
        for (Particle q : particles) {
            if (p != q) {
                double t = p.timeToHit(q);
                if (t > 1e-12) {
                    pq.add(new Event(simTime + t, p, q));
                }
            }
        }
    }

    public void simulate(double maxTime, int nEventsPerPrint) throws IOException {

        for (Particle p : particles) {
            predict(p);
        }

        pq.add(new Event(0, null, null));

        String simPath = "../outputs/sim_circular";
        File folder = new File(simPath);
        if (!folder.exists())
            folder.mkdirs();

        try (
                FileWriter outputWriter = new FileWriter(simPath + "/output.txt");) {

            outputWriter.write("# t\n");
            outputWriter.write("# x y vx vy\n");

            int countEvents = 0;

            while (!pq.isEmpty()) {
                Event event = pq.poll();
                if (!event.isValid())
                    continue;
                if (event.time < simTime - 1e-12)
                    continue;

                // Limpiar eventos inválidos periódicamente
                eventsSinceCleanup++;
                if (eventsSinceCleanup > CLEANUP_THRESHOLD && pq.size() > CLEANUP_THRESHOLD) {
                    cleanupQueue();
                }

                // mover todas las partículas
                for (Particle p : particles) {
                    p.move(event.time - simTime);
                }
                simTime = event.time;

                countEvents++;

                if (countEvents % nEventsPerPrint == 0) {
                    outputWriter.write(String.format(Locale.US, "%.6f\n", simTime));
                    for (Particle p : particles) {
                        outputWriter.write(String.format(
                                Locale.US,
                                "%.6f %.6f %.6f %.6f\n",
                                p.x, p.y, p.vx, p.vy));
                    }
                }

                if (simTime > maxTime)
                    break;

                Particle a = event.a;
                Particle b = event.b;

                // 🔥 resolver evento
                if (a != null && b != null) {

                    if (a == b) {
                        // círculo interno
                        a.bounceOffCircle(CX, CY);
                    } else {
                        // partícula-partícula
                        a.bounceOff(b);
                    }

                } else if (a != null) {
                    // pared externa
                    a.bounceOffCircle(CX, CY);
                }

                // recalcular eventos
                predict(a);
                if (b != null && b != a)
                    predict(b);
            }
        }

        System.out.println("Simulation finished.");
    }
}