package com.github.micycle1.wevo4j;

import java.util.*;
import org.locationtech.jts.algorithm.CGAlgorithmsDD;
import org.locationtech.jts.algorithm.LineIntersector;
import org.locationtech.jts.algorithm.RobustLineIntersector;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;

public final class MwvdGeom {
    public final double eps;

    public MwvdGeom(double eps) {
        if (!(eps > 0.0)) throw new IllegalArgumentException("eps must be > 0");
        this.eps = eps;
    }

    public record Circle(Coordinate c, double r) {}
    /** CCW arc from a0 to a1 on circle. */
    public record Arc(Circle circle, double a0, double a1) {}
    public record Segment(Coordinate a, Coordinate b) {}

    public double tol(double a, double b) {
        return eps * Math.max(1.0, Math.max(Math.abs(a), Math.abs(b)));
    }

    public int cmp(double a, double b) {
        double t = tol(a, b);
        if (a < b - t) return -1;
        if (a > b + t) return 1;
        return 0;
    }

    public boolean eq(double a, double b) { return cmp(a, b) == 0; }
    public boolean lt(double a, double b) { return cmp(a, b) < 0; }
    public boolean le(double a, double b) { return cmp(a, b) <= 0; }

    public boolean samePoint(Coordinate a, Coordinate b) {
        return eq(a.x, b.x) && eq(a.y, b.y);
    }

    public double dist2(Coordinate a, Coordinate b) {
        double dx = a.x - b.x, dy = a.y - b.y;
        return dx * dx + dy * dy;
    }

    public Coordinate copy(Coordinate c) { return new Coordinate(c.x, c.y); }

    public Coordinate pointOnCircle(Circle c, double a) {
        return new Coordinate(c.c.x + c.r * Math.cos(a), c.c.y + c.r * Math.sin(a));
    }

    public double angle(Coordinate center, Coordinate p) {
        double a = Math.atan2(p.y - center.y, p.x - center.x);
        return normAngle(a);
    }

    public double normAngle(double a) {
        double t = a % (2.0 * Math.PI);
        return t < 0 ? t + 2.0 * Math.PI : t;
    }

    public double ccwSweep(double a0, double a1) {
        double s = normAngle(a1) - normAngle(a0);
        return s < 0 ? s + 2.0 * Math.PI : s;
    }

    public boolean angleOnArc(double a, Arc arc) {
        double s0 = normAngle(arc.a0());
        double s1 = normAngle(arc.a1());
        double x = normAngle(a);
        double total = ccwSweep(s0, s1);
        double part = ccwSweep(s0, x);
        return le(part, total);
    }

    public boolean pointOnArc(Coordinate p, Arc arc) {
        if (!eq(dist2(p, arc.circle().c()), arc.circle().r() * arc.circle().r())) return false;
        return angleOnArc(angle(arc.circle().c(), p), arc);
    }

    public boolean pointOnSegment(Coordinate p, Segment s) {
        double dx = s.b().x - s.a().x;
        double dy = s.b().y - s.a().y;
        double len2 = dx * dx + dy * dy;
        if (eq(len2, 0.0)) return samePoint(p, s.a());
        double cross = (p.x - s.a().x) * dy - (p.y - s.a().y) * dx;
        if (Math.abs(cross) > eps * Math.max(1.0, Math.sqrt(len2))) return false;
        double dot = (p.x - s.a().x) * dx + (p.y - s.a().y) * dy;
        return le(0.0, dot) && le(dot, len2);
    }

    public int orientation(Coordinate a, Coordinate b, Coordinate c) {
        return CGAlgorithmsDD.orientationIndex(a, b, c);
    }

    public List<Coordinate> lineCircleIntersections(Coordinate o, Coordinate d, Circle c) {
        List<Coordinate> out = new ArrayList<>(2);
        double dx = d.x, dy = d.y;
        double fx = o.x - c.c().x, fy = o.y - c.c().y;
        double A = dx * dx + dy * dy;
        double B = 2.0 * (fx * dx + fy * dy);
        double C = fx * fx + fy * fy - c.r() * c.r();
        double disc = B * B - 4.0 * A * C;
        if (cmp(disc, 0.0) < 0) return out;
        if (eq(disc, 0.0)) {
            double t = -B / (2.0 * A);
            out.add(new Coordinate(o.x + t * dx, o.y + t * dy));
            return out;
        }
        double s = Math.sqrt(Math.max(0.0, disc));
        double t1 = (-B - s) / (2.0 * A);
        double t2 = (-B + s) / (2.0 * A);
        out.add(new Coordinate(o.x + t1 * dx, o.y + t1 * dy));
        addUnique(out, new Coordinate(o.x + t2 * dx, o.y + t2 * dy));
        return out;
    }

    public List<Coordinate> segmentCircleIntersections(Segment seg, Circle c) {
        Coordinate d = new Coordinate(seg.b().x - seg.a().x, seg.b().y - seg.a().y);
        List<Coordinate> raw = lineCircleIntersections(seg.a(), d, c);
        List<Coordinate> out = new ArrayList<>(2);
        for (Coordinate p : raw) if (pointOnSegment(p, seg)) addUnique(out, p);
        return out;
    }

    public List<Coordinate> circleCircleIntersections(Circle c1, Circle c2) {
        List<Coordinate> out = new ArrayList<>(2);
        double dx = c2.c().x - c1.c().x, dy = c2.c().y - c1.c().y;
        double d2 = dx * dx + dy * dy;
        double d = Math.sqrt(d2);
        if (cmp(d, 0.0) == 0) return out;
        if (cmp(d, c1.r() + c2.r()) > 0) return out;
        if (cmp(d, Math.abs(c1.r() - c2.r())) < 0) return out;

        double a = (c1.r() * c1.r() - c2.r() * c2.r() + d2) / (2.0 * d);
        double h2 = c1.r() * c1.r() - a * a;
        if (cmp(h2, 0.0) < 0) return out;
        double xm = c1.c().x + a * dx / d;
        double ym = c1.c().y + a * dy / d;

        if (eq(h2, 0.0)) {
            out.add(new Coordinate(xm, ym));
            return out;
        }

        double h = Math.sqrt(Math.max(0.0, h2));
        double rx = -dy * h / d;
        double ry =  dx * h / d;
        out.add(new Coordinate(xm + rx, ym + ry));
        addUnique(out, new Coordinate(xm - rx, ym - ry));
        return out;
    }

    public List<Coordinate> segmentSegmentIntersections(Segment s1, Segment s2) {
        List<Coordinate> out = new ArrayList<>(2);
        LineIntersector li = new RobustLineIntersector();
        li.computeIntersection(s1.a(), s1.b(), s2.a(), s2.b());
        if (!li.hasIntersection()) return out;
        for (int i = 0; i < li.getIntersectionNum(); i++) addUnique(out, li.getIntersection(i));
        if (out.isEmpty()) {
            if (pointOnSegment(s1.a(), s2)) addUnique(out, s1.a());
            if (pointOnSegment(s1.b(), s2)) addUnique(out, s1.b());
            if (pointOnSegment(s2.a(), s1)) addUnique(out, s2.a());
            if (pointOnSegment(s2.b(), s1)) addUnique(out, s2.b());
        }
        return out;
    }

    public List<Coordinate> segmentArcIntersections(Segment seg, Arc arc) {
        List<Coordinate> out = new ArrayList<>(2);
        for (Coordinate p : segmentCircleIntersections(seg, arc.circle())) {
            if (pointOnArc(p, arc)) addUnique(out, p);
        }
        return out;
    }

    public List<Coordinate> arcArcIntersections(Arc a1, Arc a2) {
        List<Coordinate> out = new ArrayList<>(2);
        for (Coordinate p : circleCircleIntersections(a1.circle(), a2.circle())) {
            if (pointOnArc(p, a1) && pointOnArc(p, a2)) addUnique(out, p);
        }
        return out;
    }

    public Coordinate rayBoxHit(Coordinate o, Coordinate dir, Envelope box) {
        double best = Double.POSITIVE_INFINITY;
        Coordinate hit = null;

        if (!eq(dir.x, 0.0)) {
            double t = (box.getMinX() - o.x) / dir.x;
            if (t > 0.0) {
                double y = o.y + t * dir.y;
                if (le(box.getMinY(), y) && le(y, box.getMaxY()) && t < best) {
                    best = t; hit = new Coordinate(box.getMinX(), y);
                }
            }
            t = (box.getMaxX() - o.x) / dir.x;
            if (t > 0.0) {
                double y = o.y + t * dir.y;
                if (le(box.getMinY(), y) && le(y, box.getMaxY()) && t < best) {
                    best = t; hit = new Coordinate(box.getMaxX(), y);
                }
            }
        }

        if (!eq(dir.y, 0.0)) {
            double t = (box.getMinY() - o.y) / dir.y;
            if (t > 0.0) {
                double x = o.x + t * dir.x;
                if (le(box.getMinX(), x) && le(x, box.getMaxX()) && t < best) {
                    best = t; hit = new Coordinate(x, box.getMinY());
                }
            }
            t = (box.getMaxY() - o.y) / dir.y;
            if (t > 0.0) {
                double x = o.x + t * dir.x;
                if (le(box.getMinX(), x) && le(x, box.getMaxX()) && t < best) {
                    best = t; hit = new Coordinate(x, box.getMaxY());
                }
            }
        }
        return hit;
    }

    private void addUnique(List<Coordinate> out, Coordinate p) {
        for (Coordinate q : out) if (samePoint(q, p)) return;
        out.add(new Coordinate(p.x, p.y));
    }
}