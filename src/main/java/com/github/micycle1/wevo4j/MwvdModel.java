package com.github.micycle1.wevo4j;

import java.util.*;
import org.locationtech.jts.geom.Coordinate;

public final class MwvdModel {
    private MwvdModel() {}

    public enum PointType { NONE, TIME, TRANS, COLL, DOM }

    public record SitePair(int a, int b) {
        public SitePair {
            if (a > b) { int t = a; a = b; b = t; }
        }
    }

    public record MovIsectId(int s1, int s2, boolean isLeft, boolean isFirst) {}
    public record ArcId(MovIsectId from, MovIsectId to) {}
    public record SwitchMark(double sqTime, Coordinate p, boolean wf) {
        public SwitchMark {
            p = new Coordinate(p.x, p.y);
        }
    }
    public record IsectPair(MovingIntersection first, MovingIntersection second) {}
    public record ExpandResult(boolean checkBothSides, boolean leftFlag) {}

    public static final class Site implements Comparable<Site> {
        public final int id;
        public final Coordinate p;
        public final double w;

        public Site(int id, Coordinate p, double w) {
            if (!(w > 0.0) || !Double.isFinite(w)) throw new IllegalArgumentException("weight must be > 0");
            this.id = id;
            this.p = new Coordinate(p.x, p.y);
            this.w = w;
        }

        public double sqTimeAt(Coordinate q) {
            double dx = q.x - p.x, dy = q.y - p.y;
            return (dx * dx + dy * dy) / (w * w);
        }

        public double angle(Coordinate q) {
            double a = Math.atan2(q.y - p.y, q.x - p.x);
            return a < 0 ? a + 2.0 * Math.PI : a;
        }

        @Override public int compareTo(Site o) {
            int c = Double.compare(w, o.w);
            return c != 0 ? c : Integer.compare(id, o.id);
        }
    }

    public static class TimePoint {
        public final Coordinate p;
        public final double sqTime;
        public TimePoint(Coordinate p, double sqTime) {
            this.p = new Coordinate(p.x, p.y);
            this.sqTime = sqTime;
        }
    }

    public static final class TransitPoint extends TimePoint {
        public final SitePair fromId;
        public final SitePair toId;
        public final PointType type;
        public TransitPoint(Coordinate p, double sqTime, SitePair fromId, SitePair toId, PointType type) {
            super(p, sqTime);
            this.fromId = fromId;
            this.toId = toId;
            this.type = type;
        }
    }

    public static abstract class TrajectorySection {
        protected final MwvdGeom g;
        protected final TransitPoint start;
        protected final TransitPoint end;
        protected final Site site1;
        protected final Site site2;

        protected TrajectorySection(MwvdGeom g, TransitPoint start, TransitPoint end, Site site1, Site site2) {
            this.g = g; this.start = start; this.end = end; this.site1 = site1; this.site2 = site2;
        }

        public boolean includes(double t) { return g.le(start.sqTime, t) && g.le(t, end.sqTime); }
        public abstract Coordinate pointAt(double t);
        public abstract boolean supportContains(Coordinate p);
        public abstract List<TimePoint> intersections(TrajectorySection other);
    }

    public static final class PointPointSection extends TrajectorySection {
        public final boolean equal;
        public final MwvdGeom.Arc arc;
        public final MwvdGeom.Segment seg;

        public PointPointSection(MwvdGeom g, TransitPoint start, TransitPoint end, Site s1, Site s2, MwvdGeom.Arc arc) {
            super(g, start, end, s1, s2);
            this.equal = false; this.arc = arc; this.seg = null;
        }

        public PointPointSection(MwvdGeom g, TransitPoint start, TransitPoint end, Site s1, Site s2, MwvdGeom.Segment seg) {
            super(g, start, end, s1, s2);
            this.equal = true; this.arc = null; this.seg = seg;
        }

        @Override public Coordinate pointAt(double t) {
            if (g.eq(t, start.sqTime)) return new Coordinate(start.p.x, start.p.y);
            if (g.eq(t, end.sqTime)) return new Coordinate(end.p.x, end.p.y);
            MwvdGeom.Circle grow = new MwvdGeom.Circle(site1.p, site1.w * Math.sqrt(Math.max(0.0, t)));
            List<Coordinate> cand = equal ? g.segmentCircleIntersections(seg, grow)
                                          : g.circleCircleIntersections(arc.circle(), grow);
            for (Coordinate p : cand) {
                if ((equal ? g.pointOnSegment(p, seg) : g.pointOnArc(p, arc)) && g.eq(site1.sqTimeAt(p), t)) {
                    return p;
                }
            }
            Coordinate best = null;
            double bestErr = Double.POSITIVE_INFINITY;
            for (Coordinate p : cand) {
                double err = Math.abs(site1.sqTimeAt(p) - t);
                if (err < bestErr) { bestErr = err; best = p; }
            }
            assert best != null;
            return best;
        }

        @Override public boolean supportContains(Coordinate p) {
            return equal ? g.pointOnSegment(p, seg) : g.pointOnArc(p, arc);
        }

        @Override public List<TimePoint> intersections(TrajectorySection other) {
            List<TimePoint> out = new ArrayList<>();
            PointPointSection o = (PointPointSection) other;
            List<Coordinate> cand;
            if (equal && o.equal) cand = g.segmentSegmentIntersections(seg, o.seg);
            else if (equal) cand = g.segmentArcIntersections(seg, o.arc);
            else if (o.equal) cand = g.segmentArcIntersections(o.seg, arc);
            else cand = g.arcArcIntersections(arc, o.arc);

            for (Coordinate p : cand) {
                double t = site1.sqTimeAt(p);
                if (includes(t) && o.includes(t)) addUnique(out, new TimePoint(p, t), g);
            }
            return out;
        }

        private static void addUnique(List<TimePoint> out, TimePoint tp, MwvdGeom g) {
            for (TimePoint q : out) {
                if (g.eq(q.sqTime, tp.sqTime) && g.samePoint(q.p, tp.p)) return;
            }
            out.add(tp);
        }
    }

    public static final class Trajectory {
        public final Site site1;
        public final Site site2;
        public final boolean isLeft;
        public final boolean isFirst;
        public final List<TrajectorySection> secs = new ArrayList<>();

        public Trajectory(Site site1, Site site2, boolean isLeft, boolean isFirst) {
            this.site1 = site1; this.site2 = site2; this.isLeft = isLeft; this.isFirst = isFirst;
        }

        public SitePair id() { return new SitePair(site1.id, site2.id); }
        public void add(TrajectorySection s) { secs.add(s); }
        public TransitPoint start() { return (TransitPoint) secs.get(0).start; }
        public TransitPoint end() { return (TransitPoint) secs.get(secs.size() - 1).end; }

        public Coordinate pointAt(double t) {
            for (TrajectorySection s : secs) if (s.includes(t)) return s.pointAt(t);
            return secs.get(secs.size() - 1).pointAt(end().sqTime);
        }

        public List<TimePoint> intersections(Trajectory other) {
            if (id().equals(other.id())) return List.of();
            List<TimePoint> out = new ArrayList<>();
            for (TrajectorySection a : secs) for (TrajectorySection b : other.secs) {
                for (TimePoint tp : a.intersections(b)) addUnique(out, tp, ((PointPointSection) a).g);
            }
            return out;
        }

        public Site otherSite(int siteId) {
            assert site1.id == siteId || site2.id == siteId;
            return site1.id == siteId ? site2 : site1;
        }

        private static void addUnique(List<TimePoint> out, TimePoint tp, MwvdGeom g) {
            for (TimePoint q : out) {
                if (g.eq(q.sqTime, tp.sqTime) && g.samePoint(q.p, tp.p)) return;
            }
            out.add(tp);
        }
    }

    public static final class Bisector {
        public final MwvdGeom g;
        public final Site site1;   // always >= site2 by weight/id
        public final Site site2;
        public final Coordinate lineA; // site2 -> site1
        public final Coordinate lineB;
        public final List<Trajectory> trajs = new ArrayList<>(2);

        public Bisector(MwvdGeom g, Site a, Site b, double tCap) {
            this.g = g;
            if (a.compareTo(b) < 0) { Site t = a; a = b; b = t; }
            this.site1 = a; this.site2 = b;
            this.lineA = site2.p;
            this.lineB = site1.p;
            compute(tCap);
        }

        public SitePair id() { return new SitePair(site1.id, site2.id); }

        private void compute(double tCap) {
            if (g.eq(site1.w, site2.w)) {
                Coordinate mid = new Coordinate((site1.p.x + site2.p.x) * 0.5, (site1.p.y + site2.p.y) * 0.5);
                double dx = site1.p.x - site2.p.x, dy = site1.p.y - site2.p.y;
                double n = Math.hypot(dx, dy);
                Coordinate dLeft = new Coordinate(-dy / n, dx / n);
                Coordinate dRight = new Coordinate(-dLeft.x, -dLeft.y);

                double t0 = site1.sqTimeAt(mid);
                double rr = site1.w * site1.w * tCap;
                double base = (mid.x - site1.p.x) * (mid.x - site1.p.x) + (mid.y - site1.p.y) * (mid.y - site1.p.y);
                double s = Math.sqrt(Math.max(0.0, rr - base));

                Coordinate ePos = new Coordinate(mid.x + dLeft.x * s, mid.y + dLeft.y * s);
                Coordinate eNeg = new Coordinate(mid.x + dRight.x * s, mid.y + dRight.y * s);

                TransitPoint coll = new TransitPoint(mid, t0, id(), id(), PointType.COLL);
                TransitPoint domPos = new TransitPoint(ePos, tCap, id(), id(), PointType.DOM);
                TransitPoint domNeg = new TransitPoint(eNeg, tCap, id(), id(), PointType.DOM);

                Trajectory pos = new Trajectory(site1, site2, true, true);
                Trajectory neg = new Trajectory(site1, site2, false, true);
                pos.add(new PointPointSection(g, coll, domPos, site1, site2, new MwvdGeom.Segment(mid, ePos)));
                neg.add(new PointPointSection(g, coll, domNeg, site1, site2, new MwvdGeom.Segment(mid, eNeg)));
                trajs.add(pos); trajs.add(neg);
                return;
            }

            double w1s = site1.w * site1.w, w2s = site2.w * site2.w;
            double d2 = g.dist2(site1.p, site2.p);
            double denom = w1s - w2s;
            double cx = (w1s * site2.p.x - w2s * site1.p.x) / denom;
            double cy = (w1s * site2.p.y - w2s * site1.p.y) / denom;
            double r2 = (w1s * w2s * d2) / ((w2s - w1s) * (w2s - w1s));
            MwvdGeom.Circle circ = new MwvdGeom.Circle(new Coordinate(cx, cy), Math.sqrt(Math.max(0.0, r2)));

            List<Coordinate> cuts = g.lineCircleIntersections(site2.p,
                    new Coordinate(site1.p.x - site2.p.x, site1.p.y - site2.p.y), circ);
            assert cuts.size() == 2;
            if (cuts.size() != 2) return;

            Coordinate p1 = cuts.get(0), p2 = cuts.get(1);
            double t1 = site1.sqTimeAt(p1), t2 = site1.sqTimeAt(p2);
            if (g.lt(t2, t1)) {
                Coordinate tp = p1; p1 = p2; p2 = tp;
                double tt = t1; t1 = t2; t2 = tt;
            }

            TransitPoint coll = new TransitPoint(p1, t1, id(), id(), PointType.COLL);
            TransitPoint dom  = new TransitPoint(p2, t2, id(), id(), PointType.DOM);

            MwvdGeom.Arc a = new MwvdGeom.Arc(circ, g.angle(circ.c(), p1), g.angle(circ.c(), p2));
            MwvdGeom.Arc b = new MwvdGeom.Arc(circ, g.angle(circ.c(), p2), g.angle(circ.c(), p1));

            Coordinate sa = sample(a), sb = sample(b);
            boolean aPos = g.orientation(lineA, lineB, sa) > 0;
            MwvdGeom.Arc posArc = aPos ? a : b;
            MwvdGeom.Arc negArc = aPos ? b : a;

            Trajectory pos = new Trajectory(site1, site2, true, true);
            Trajectory neg = new Trajectory(site1, site2, false, true);
            pos.add(new PointPointSection(g, coll, dom, site1, site2, posArc));
            neg.add(new PointPointSection(g, coll, dom, site1, site2, negArc));
            trajs.add(pos); trajs.add(neg);
        }

        private Coordinate sample(MwvdGeom.Arc a) {
            double m = g.normAngle(a.a0() + 0.5 * g.ccwSweep(a.a0(), a.a1()));
            return g.pointOnCircle(a.circle(), m);
        }

        public Trajectory findTraj(Coordinate p) {
            int side = g.orientation(lineA, lineB, p);
            if (side > 0) return trajs.get(0);
            if (side < 0) return trajs.get(1);
            for (Trajectory t : trajs) {
                if (((PointPointSection) t.secs.get(0)).supportContains(p)) return t;
            }
            assert false;
            return trajs.get(0);
        }
    }

    public static final class MovingIntersection {
        public final Trajectory traj;
        public boolean wf = true;
        public final List<SwitchMark> switches = new ArrayList<>();

        public MovingIntersection(Trajectory traj) { this.traj = traj; }
        public MovIsectId id() {
            SitePair id = traj.id();
            return new MovIsectId(id.a(), id.b(), traj.isLeft, traj.isFirst);
        }
        public Coordinate pointAt(double t) { return traj.pointAt(t); }
        
        public void setIsWfVert(double t, Coordinate p, boolean value) {
            wf = value;
            switches.add(new SwitchMark(t, p, value));
        }
    }

    public static final class OffCircle {
        public final MwvdGeom g;
        public final Site site;
        public boolean active = true;
        public final Map<ArcId, Boolean> arcs = new HashMap<>();
        public final Map<MovIsectId, MovingIntersection> isects = new HashMap<>();
        public final Map<MovIsectId, MovIsectId> lefts = new HashMap<>();
        public final Map<MovIsectId, MovIsectId> rights = new HashMap<>();

        public OffCircle(MwvdGeom g, Site site) { this.g = g; this.site = site; }

        public boolean includesIsect(MovingIntersection i) { return isects.containsKey(i.id()); }
        public boolean isActive() { return active; }

        public MovingIntersection neighbor(MovingIntersection i, boolean left) {
            MovIsectId id = left ? lefts.get(i.id()) : rights.get(i.id());
            return id == null ? null : isects.get(id);
        }

        public MovingIntersection searchNeighbor(MovingIntersection i, boolean left) { return neighbor(i, left); }

        public MovingIntersection searchNeighbor(double t, Coordinate p, boolean left) {
            double a = site.angle(p);
            MovingIntersection res = searchNeighborByAngle(t, a, left);
            if (res == null) res = searchNeighborByAngle(t, left ? 2.0 * Math.PI : 0.0, left);
            return res;
        }

        private MovingIntersection searchNeighborByAngle(double t, double angle, boolean left) {
            double best = Double.POSITIVE_INFINITY;
            MovingIntersection out = null;
            for (MovingIntersection i : isects.values()) {
                double ai = site.angle(i.pointAt(t));
                boolean ok = left ? ai < angle : ai > angle;
                double diff = Math.abs(angle - ai);
                if (ok && diff < best) { best = diff; out = i; }
            }
            return out;
        }

        public boolean isInActiveArc(double t, Coordinate p) {
            if (isects.isEmpty()) return active;
            MovingIntersection left = searchNeighbor(t, p, true);
            MovingIntersection right = searchNeighbor(t, p, false);
            if (left == null || right == null) return false;
            return containsArc(new ArcId(left.id(), right.id()), false);
        }

        public void spawnArc(double t, Coordinate p, MovingIntersection i1, MovingIntersection i2, boolean isActive, boolean pierces) {
            if (!active) return;
            MovingIntersection from = isActive ? i1 : i2;
            MovingIntersection to = isActive ? i2 : i1;
            if (!isActive && pierces) { MovingIntersection tmp = from; from = to; to = tmp; }

            if (isects.isEmpty()) {
                if (isActive) insertArc(from.id(), to.id(), false);
                insertArc(to.id(), from.id(), true);
                insertIsect(i1); i1.setIsWfVert(t, p, true);
                insertIsect(i2); i2.setIsWfVert(t, p, true);
                return;
            }

            MovingIntersection left = searchNeighbor(t, p, true);
            MovingIntersection right = searchNeighbor(t, p, false);
            assert left != null && right != null;

            ArcId old = new ArcId(left.id(), right.id());
            if (containsArc(old, true)) {
                boolean onWf = arcs.get(old);
                eraseArc(old, true);
                insertArc(left.id(), from.id(), onWf);
                if (isActive) insertArc(from.id(), to.id(), pierces);
                insertArc(to.id(), right.id(), onWf);
                insertIsect(i1); i1.setIsWfVert(t, p, onWf || (isActive && pierces));
                insertIsect(i2); i2.setIsWfVert(t, p, onWf || (isActive && pierces));
            }
        }

        public IsectPair deleteArc(double t, Coordinate p, MovingIntersection from, MovingIntersection to, boolean isActive) {
            IsectPair newArc = null;

            if (arcs.size() > 2) {
                MovingIntersection left = searchNeighbor(from, true);
                MovingIntersection right = searchNeighbor(to, false);

                if (left != null && right != null) {
                    ArcId a = new ArcId(left.id(), from.id());
                    ArcId b = new ArcId(from.id(), to.id());
                    ArcId c = new ArcId(to.id(), right.id());

                    if (containsArc(a, false) && containsArc(b, false) && containsArc(c, false)) {
                        boolean onWf = arcs.get(a);
                        eraseArc(a, false);
                        eraseArc(b, false);
                        eraseArc(c, false);
                        insertArc(left.id(), right.id(), onWf);
                        newArc = new IsectPair(left, right);
                    } else if (containsArc(new ArcId(from.id(), to.id()), false)) {
                        eraseArc(new ArcId(from.id(), to.id()), false);
                    } else if (containsArc(new ArcId(to.id(), from.id()), false)) {
                        eraseArc(new ArcId(to.id(), from.id()), false);
                    }
                } else if (containsArc(new ArcId(from.id(), to.id()), false)) {
                    eraseArc(new ArcId(from.id(), to.id()), false);
                } else if (containsArc(new ArcId(to.id(), from.id()), false)) {
                    eraseArc(new ArcId(to.id(), from.id()), false);
                }
            } else {
                eraseArc(new ArcId(from.id(), to.id()), false);
                eraseArc(new ArcId(to.id(), from.id()), false);
                if (arcs.isEmpty()) active = isActive;
            }

            if (isects.remove(from.id()) != null) from.setIsWfVert(t, p, false);
            if (isects.remove(to.id()) != null) to.setIsWfVert(t, p, false);

            return newArc;
        }

        public IsectPair deleteArcUnordered(double t, Coordinate p, MovingIntersection a, MovingIntersection b, boolean isActive) {
            if (containsArc(new ArcId(a.id(), b.id()), false)) return deleteArc(t, p, a, b, isActive);
            if (containsArc(new ArcId(b.id(), a.id()), false)) return deleteArc(t, p, b, a, isActive);
            return null;
        }

        private ArcId findIncidentArc(MovIsectId id) {
            MovIsectId left = lefts.get(id);
            if (left != null && containsArc(new ArcId(left, id), false)) return new ArcId(left, id);
            MovIsectId right = rights.get(id);
            if (right != null && containsArc(new ArcId(id, right), false)) return new ArcId(id, right);
            return null;
        }


        public boolean collapseArc(double t, MovingIntersection from, MovingIntersection to) {
            boolean less = containsArc(new ArcId(from.id(), to.id()), false);
            MovingIntersection n = searchNeighbor(from, less);
            assert n != null;

            ArcId old1 = less ? new ArcId(from.id(), to.id()) : new ArcId(to.id(), from.id());
            ArcId old2 = less ? new ArcId(n.id(), from.id()) : new ArcId(from.id(), n.id());
            ArcId neu  = less ? new ArcId(n.id(), to.id()) : new ArcId(to.id(), n.id());

            eraseArc(old1, true);
            if (containsArc(old2, false)) {
                boolean onWf = arcs.get(old2);
                eraseArc(old2, true);
                isects.remove(from.id());
                insertArc(neu.from(), neu.to(), onWf);
            } else {
                isects.remove(from.id());
                isects.remove(to.id());
            }
            return less;
        }

        public ExpandResult expandIsect(double t, Coordinate evPnt, MovingIntersection from, MovingIntersection to, boolean onWf1) {
            int side = g.orientation(site.p, evPnt, from.traj.otherSite(site.id).p);
            boolean less = side < 0;
            MovingIntersection n = searchNeighbor(from, !less);
            ArcId new1 = less ? new ArcId(from.id(), to.id()) : new ArcId(to.id(), from.id());

            if (n != null) {
                ArcId old = less ? new ArcId(from.id(), n.id()) : new ArcId(n.id(), from.id());
                ArcId new2 = less ? new ArcId(to.id(), n.id()) : new ArcId(n.id(), to.id());
                if (containsArc(old, true)) {
                    boolean onWf2 = arcs.get(old);
                    eraseArc(old, true);
                    insertArc(new1.from(), new1.to(), onWf1);
                    insertArc(new2.from(), new2.to(), onWf2);
                }
            } else {
                insertArc(new1.from(), new1.to(), onWf1);
            }

            insertIsect(to);
            return new ExpandResult(n != null, less);
        }

        public boolean replaceIsect(double t, MovingIntersection oldIsect, MovingIntersection newIsect) {
            ArcId arcId = findIncidentArc(oldIsect.id());
            if (arcId == null) return false;

            boolean less = oldIsect.id().equals(arcId.from());
            MovIsectId other = less ? arcId.to() : arcId.from();
            if (!isects.containsKey(other)) return false;

            boolean onWf = oldIsect.wf;
            isects.remove(oldIsect.id());
            insertIsect(newIsect);

            if (containsArc(arcId, false)) {
                eraseArc(arcId, false);
                insertArc(less ? newIsect.id() : other, less ? other : newIsect.id(), onWf);
            }
            return !less;
        }

        private void insertArc(MovIsectId from, MovIsectId to, boolean onWf) {
            arcs.put(new ArcId(from, to), onWf);
            rights.put(from, to);
            lefts.put(to, from);
        }

        private void eraseArc(ArcId id, boolean assertPresent) {
            if (assertPresent) assert arcs.containsKey(id);
            arcs.remove(id);
            rights.remove(id.from());
            lefts.remove(id.to());
        }

        private void insertIsect(MovingIntersection i) { isects.put(i.id(), i); }
        private boolean containsArc(ArcId id, boolean assertPresent) {
            if (assertPresent) assert arcs.containsKey(id);
            return arcs.containsKey(id);
        }
    }
}
