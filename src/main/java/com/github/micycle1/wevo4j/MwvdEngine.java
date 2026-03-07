package com.github.micycle1.wevo4j;

import java.util.*;
import org.locationtech.jts.geom.*;

public final class MwvdEngine {
    private final MwvdGeom g;
    private final GeometryFactory gf;
    private final Envelope clip;
    private final double chordTol;

    private final PriorityQueue<MwvdEvents.Event> queue;
    private final Map<MwvdModel.SitePair, MwvdModel.Bisector> bisecs = new HashMap<>();
    private final Map<Integer, MwvdModel.OffCircle> offCircs = new HashMap<>();
    private final Map<MwvdModel.MovIsectId, MwvdModel.MovingIntersection> isects = new HashMap<>();

    public MwvdEngine(double eps) {
        this(eps, new GeometryFactory(), null, 1e-2);
    }

    public MwvdEngine(double eps, GeometryFactory gf, Envelope clip, double chordTol) {
        this.g = new MwvdGeom(eps);
        this.gf = gf;
        this.clip = clip;
        this.chordTol = chordTol;
        this.queue = MwvdEvents.newQueue(g);
    }

    public GeometryCollection computeRegions(List<Coordinate> input) {
        List<MwvdModel.Site> sites = buildSites(input);
        if (sites.isEmpty()) return gf.createGeometryCollection();

        Envelope box = this.clip != null ? new Envelope(this.clip) : autoClip(sites);
        double tCap = capTime(sites, box);

        for (MwvdModel.Site s : sites) {
            offCircs.put(s.id, new MwvdModel.OffCircle(g, s));
        }

        for (int i = 0; i < sites.size(); i++) {
            for (int j = i + 1; j < sites.size(); j++) {
                MwvdModel.Bisector b = bisector(sites.get(i), sites.get(j), tCap);
                if (b.trajs.size() != 2) continue;

                MwvdModel.Trajectory t1 = b.trajs.get(0);
                MwvdModel.Trajectory t2 = b.trajs.get(1);

                queue.add(new MwvdEvents.CollEvent(t1.start().p, t1.start().sqTime, sites.get(i), t1, t2, false));
                queue.add(new MwvdEvents.CollEvent(t1.start().p, t1.start().sqTime, sites.get(j), t1, t2, false));
            }
        }

        run();

        List<MwvdExporter.EdgePiece> edges = collectEdges();
        return new MwvdExporter(g, gf, chordTol).exportRegions(edges, sites, box);
    }

    private List<MwvdModel.Site> buildSites(List<Coordinate> input) {
        List<MwvdModel.Site> sites = new ArrayList<>(input.size());
        for (int i = 0; i < input.size(); i++) {
            Coordinate c = input.get(i);
            double w = Double.isFinite(c.getZ()) ? c.getZ() : 1.0;
            sites.add(new MwvdModel.Site(i, c, w));
        }
        sites.sort((a, b) -> -a.compareTo(b));
        return sites;
    }

    private Envelope autoClip(List<MwvdModel.Site> sites) {
        Envelope e = new Envelope();
        for (MwvdModel.Site s : sites) e.expandToInclude(s.p);
        double span = Math.max(e.getWidth(), e.getHeight());
        if (span == 0.0) span = 1.0;
        e.expandBy(8.0 * span);
        return e;
    }

    private double capTime(List<MwvdModel.Site> sites, Envelope box) {
        Coordinate[] cs = new Coordinate[]{
                new Coordinate(box.getMinX(), box.getMinY()),
                new Coordinate(box.getMaxX(), box.getMinY()),
                new Coordinate(box.getMaxX(), box.getMaxY()),
                new Coordinate(box.getMinX(), box.getMaxY())
        };
        double max = 1.0;
        for (MwvdModel.Site s : sites) for (Coordinate c : cs) max = Math.max(max, s.sqTimeAt(c));
        return 16.0 * max;
    }

    private MwvdModel.Bisector bisector(MwvdModel.Site a, MwvdModel.Site b, double tCap) {
        MwvdModel.SitePair id = new MwvdModel.SitePair(a.id, b.id);
        return bisecs.computeIfAbsent(id, k -> new MwvdModel.Bisector(g, a, b, tCap));
    }

    private void run() {
        while (!queue.isEmpty()) {
            MwvdEvents.Event ev = queue.poll();
            switch (ev.type()) {
                case COLL -> handle((MwvdEvents.CollEvent) ev);
                case DOM  -> handle((MwvdEvents.DomEvent) ev);
                case EDGE -> handle((MwvdEvents.EdgeEvent) ev);
            }
        }
    }

    private void handle(MwvdEvents.CollEvent ev1) {
        MwvdEvents.Event raw = queue.poll();
        assert raw instanceof MwvdEvents.CollEvent;
        MwvdEvents.CollEvent ev2 = (MwvdEvents.CollEvent) raw;

        boolean dom = ev1.site.compareTo(ev2.site) > 0;
        boolean valid = isValidColl(ev1, ev2);

        process(ev1, dom, valid, ev1.pierces);
        process(ev2, !dom, valid, ev1.pierces);
    }

    private boolean isValidColl(MwvdEvents.CollEvent a, MwvdEvents.CollEvent b) {
        return offCircs.get(a.site.id).isInActiveArc(a.sqTime, a.p)
            && offCircs.get(b.site.id).isInActiveArc(b.sqTime, b.p);
    }

    private void process(MwvdEvents.CollEvent ev, boolean dom, boolean valid, boolean pierces) {
        MwvdModel.OffCircle off = offCircs.get(ev.site.id);
        MwvdModel.MovingIntersection i1 = makeMovIsect(ev.traj1);
        MwvdModel.MovingIntersection i2 = makeMovIsect(ev.traj2);
        assert ev.traj1.isLeft && !ev.traj2.isLeft;

        if (valid) {
            off.spawnArc(ev.sqTime, i1, i2, dom, pierces);
            checkEdgeEv(off, ev.sqTime, i1, pierces ? !dom : dom);
            checkEdgeEv(off, ev.sqTime, i2, pierces ? dom : !dom);
        }

        queue.add(new MwvdEvents.DomEvent(i1.traj.end().p, i1.traj.end().sqTime, ev.site, i1, i2));
    }

    private void handle(MwvdEvents.DomEvent ev1) {
        MwvdEvents.Event raw = queue.poll();
        assert raw instanceof MwvdEvents.DomEvent;
        MwvdEvents.DomEvent ev2 = (MwvdEvents.DomEvent) raw;

        boolean dom = ev1.site.compareTo(ev2.site) > 0;
        if (isValidDom(ev1, ev2)) {
            process(ev1, dom);
            process(ev2, !dom);
        }
    }

    private boolean isValidDom(MwvdEvents.DomEvent a, MwvdEvents.DomEvent b) {
        MwvdModel.OffCircle o1 = offCircs.get(a.site.id), o2 = offCircs.get(b.site.id);
        return (o1.isInActiveArc(a.sqTime, a.p) && o2.isInActiveArc(a.sqTime, a.p))
            || (o1.isActive() && o2.isActive());
    }

    private void process(MwvdEvents.DomEvent ev, boolean dom) {
        MwvdModel.OffCircle off = offCircs.get(ev.site.id);
        MwvdModel.IsectPair pair = off.deleteArc(ev.sqTime, ev.isect1, ev.isect2, dom);
        if (pair != null) checkEdgeEv(off, ev.sqTime, pair.first(), pair.second());
    }

    private void handle(MwvdEvents.EdgeEvent edgeEv) {
        double t = edgeEv.sqTime;
        MwvdModel.Site site1 = edgeEv.site;
        MwvdModel.Site site2 = otherSite(site1, edgeEv.isect1);
        MwvdModel.Site site3 = otherSite(site1, edgeEv.isect2);

        Set<Integer> seen = new HashSet<>();
        seen.add(site1.id);
        List<MwvdEvents.EdgeEvent> twins = new ArrayList<>();

        while (!queue.isEmpty() && queue.peek() instanceof MwvdEvents.EdgeEvent nxt &&
                g.eq(edgeEv.sqTime, nxt.sqTime) && g.samePoint(edgeEv.p, nxt.p)) {
            MwvdEvents.EdgeEvent twin = (MwvdEvents.EdgeEvent) queue.poll();
            MwvdModel.OffCircle off = offCircs.get(twin.site.id);
            if (!seen.contains(twin.site.id) && off.includesIsect(twin.isect1) && off.includesIsect(twin.isect2)) {
                twins.add(twin);
                seen.add(twin.site.id);
            }
        }

        if (!twins.isEmpty()) {
            List<MwvdModel.Site> ss = new ArrayList<>(List.of(site1, site2, site3));
            ss.sort(MwvdModel.Site::compareTo);
            deleteLowestArc(edgeEv, ss.get(0), ss.get(1), ss.get(2));
            return;
        }

        if (site1.compareTo(site2) > 0 && site1.compareTo(site3) > 0) {
            MwvdModel.Site low = site2.compareTo(site3) < 0 ? site2 : site3;
            MwvdModel.Site med = site2.compareTo(site3) < 0 ? site3 : site2;

            MwvdModel.MovingIntersection i1 = makeMovIsectAt(edgeEv.p, low, med);
            MwvdModel.MovingIntersection i2 = makeMovIsectAt(edgeEv.p, low, site1);
            MwvdModel.MovingIntersection i3 = makeMovIsectAt(edgeEv.p, med, site1);
            if (i1 == null || i2 == null || i3 == null) return;

            MwvdModel.OffCircle lowC = offCircs.get(low.id), medC = offCircs.get(med.id), hiC = offCircs.get(site1.id);
            if (!lowC.includesIsect(i2) || !medC.includesIsect(i3) || !hiC.includesIsect(i2) || !hiC.includesIsect(i3)) return;

            boolean wf2 = i2.wf, wf3 = i3.wf;
            MwvdModel.ExpandResult chk = medC.expandIsect(t, edgeEv.p, i3, i1, wf2 && !wf3);
            boolean left1 = hiC.collapseArc(t, i2, i3);
            boolean left3 = lowC.replaceIsect(t, i2, i1);

            checkEdgeEv(hiC, t, i3, left1);
            if (chk.checkBothSides()) { checkEdgeEv(medC, t, i1, true); checkEdgeEv(medC, t, i1, false); }
            else checkEdgeEv(medC, t, i1, chk.leftFlag());
            checkEdgeEv(lowC, t, i1, left3);

            if (wf2 && wf3) { i1.setIsWfVert(t, true); i2.setIsWfVert(t, false); i3.setIsWfVert(t, false); }
            else if (wf2 && !wf3) { i1.setIsWfVert(t, true); i2.setIsWfVert(t, false); i3.setIsWfVert(t, true); }
            else if (!wf2 && !wf3) i1.setIsWfVert(t, false);
            return;
        }

        if (!(site1.compareTo(site2) < 0 && site1.compareTo(site3) < 0)) {
            MwvdModel.Site low = site2.compareTo(site3) < 0 ? site2 : site3;
            MwvdModel.Site high = site2.compareTo(site3) < 0 ? site3 : site2;

            MwvdModel.MovingIntersection i1 = makeMovIsectAt(edgeEv.p, low, site1);
            MwvdModel.MovingIntersection i2 = makeMovIsectAt(edgeEv.p, low, high);
            MwvdModel.MovingIntersection i3 = makeMovIsectAt(edgeEv.p, site1, high);
            if (i1 == null || i2 == null || i3 == null) return;

            MwvdModel.OffCircle lowC = offCircs.get(low.id), medC = offCircs.get(site1.id), hiC = offCircs.get(high.id);
            if (!lowC.includesIsect(i1) || !medC.includesIsect(i1) || !medC.includesIsect(i3) || !hiC.includesIsect(i3)) return;

            boolean wf1 = i1.wf, wf3 = i3.wf;
            MwvdModel.ExpandResult chk = hiC.expandIsect(t, edgeEv.p, i3, i2, wf1 && !wf3);
            boolean left2 = medC.collapseArc(t, i1, i3);
            boolean left3 = lowC.replaceIsect(t, i1, i2);

            if (chk.checkBothSides()) { checkEdgeEv(hiC, t, i2, true); checkEdgeEv(hiC, t, i2, false); }
            else checkEdgeEv(hiC, t, i2, chk.leftFlag());

            checkEdgeEv(medC, t, i3, left2);
            checkEdgeEv(lowC, t, i2, left3);

            if (wf1 && wf3) { i1.setIsWfVert(t, false); i2.setIsWfVert(t, true); i3.setIsWfVert(t, false); }
            else if (wf1 && !wf3) { i1.setIsWfVert(t, false); i2.setIsWfVert(t, true); i3.setIsWfVert(t, true); }
            else if (!wf1 && !wf3) i2.setIsWfVert(t, false);
        }
    }

    private boolean deleteLowestArc(MwvdEvents.EdgeEvent edgeEv,
                                    MwvdModel.Site lowSite, MwvdModel.Site medSite, MwvdModel.Site highSite) {
        double t = edgeEv.sqTime;
        MwvdModel.MovingIntersection i1 = makeMovIsectAt(edgeEv.p, lowSite, medSite);
        MwvdModel.MovingIntersection i2 = makeMovIsectAt(edgeEv.p, lowSite, highSite);
        MwvdModel.MovingIntersection i3 = makeMovIsectAt(edgeEv.p, medSite, highSite);
        if (i1 == null || i2 == null || i3 == null) return false;

        MwvdModel.OffCircle low = offCircs.get(lowSite.id), med = offCircs.get(medSite.id), hi = offCircs.get(highSite.id);
        if (!low.includesIsect(i1) || !low.includesIsect(i2) || !med.includesIsect(i1) || !med.includesIsect(i3)
                || !hi.includesIsect(i2) || !hi.includesIsect(i3)) return false;

        boolean wf1 = i1.wf, wf2 = i2.wf, wf3 = i3.wf;
        boolean left1 = hi.collapseArc(t, i2, i3);
        boolean left2 = med.collapseArc(t, i1, i3);
        low.deleteArcUnordered(t, i1, i2, false);

        checkEdgeEv(hi, t, i3, left1);
        checkEdgeEv(med, t, i3, left2);

        if (wf1 && wf2 && wf3) { i1.setIsWfVert(t, false); i2.setIsWfVert(t, false); i3.setIsWfVert(t, false); }
        else if (wf1 && wf2) { i1.setIsWfVert(t, false); i2.setIsWfVert(t, false); i3.setIsWfVert(t, true); }
        return true;
    }

    private MwvdModel.Site otherSite(MwvdModel.Site site, MwvdModel.MovingIntersection i) {
        return i.traj.site1.id == site.id ? i.traj.site2 : i.traj.site1;
    }

    private MwvdModel.MovingIntersection makeMovIsectAt(Coordinate p, MwvdModel.Site s1, MwvdModel.Site s2) {
        MwvdModel.Bisector b = bisecs.get(new MwvdModel.SitePair(s1.id, s2.id));
        if (b == null) return null;
        MwvdModel.Trajectory t = b.findTraj(p);
        return makeMovIsect(t);
    }

    private MwvdModel.MovingIntersection makeMovIsect(MwvdModel.Trajectory t) {
        MwvdModel.MovIsectId id = new MwvdModel.MovIsectId(t.id().a(), t.id().b(), t.isLeft, t.isFirst);
        return isects.computeIfAbsent(id, k -> new MwvdModel.MovingIntersection(t));
    }

    private void checkEdgeEv(MwvdModel.OffCircle off, double tNow, MwvdModel.MovingIntersection i, boolean left) {
        MwvdModel.MovingIntersection other = off.neighbor(i, left);
        if (other != null) checkEdgeEv(off, tNow, i, other);
    }

    private void checkEdgeEv(MwvdModel.OffCircle off, double tNow,
                             MwvdModel.MovingIntersection i, MwvdModel.MovingIntersection other) {
        MwvdModel.TimePoint best = null;
        for (MwvdModel.TimePoint tp : i.traj.intersections(other.traj)) {
            if (g.lt(tNow, tp.sqTime) && (best == null || g.lt(tp.sqTime, best.sqTime))) best = tp;
        }
        if (best != null) queue.add(new MwvdEvents.EdgeEvent(best.p, best.sqTime, off.site, i, other));
    }

    private List<MwvdExporter.EdgePiece> collectEdges() {
        List<MwvdExporter.EdgePiece> out = new ArrayList<>();

        for (MwvdModel.MovingIntersection mi : isects.values()) {
            if (mi.switches.isEmpty()) continue;

            List<MwvdModel.SwitchMark> sw = new ArrayList<>(mi.switches);
            double endT = mi.traj.end().sqTime;

            if (sw.get(sw.size() - 1).wf() && g.lt(sw.get(sw.size() - 1).sqTime(), endT)) {
                sw.add(new MwvdModel.SwitchMark(endT, false));
            }

            for (int i = 0; i + 1 < sw.size(); i++) {
                double t1 = sw.get(i).sqTime();
                double t2 = sw.get(i + 1).sqTime();

                if (!sw.get(i).wf() || !g.lt(t1, t2)) continue;

                double tm = 0.5 * (t1 + t2);
                MwvdModel.TrajectorySection sec = sectionAt(mi.traj, tm);
                if (!(sec instanceof MwvdModel.PointPointSection pp)) continue;

                Coordinate p1 = mi.pointAt(t1);
                Coordinate p2 = mi.pointAt(t2);

                if (pp.equal) {
                    out.add(new MwvdExporter.SegPiece(
                            new MwvdGeom.Segment(p1, p2),
                            pp.site1.id, pp.site2.id));
                } else {
                    double a1 = g.angle(pp.arc.circle().c(), p1);
                    double a2 = g.angle(pp.arc.circle().c(), p2);

                    MwvdGeom.Arc arc = mi.traj.isLeft
                            ? new MwvdGeom.Arc(pp.arc.circle(), a1, a2)
                            : new MwvdGeom.Arc(pp.arc.circle(), a2, a1);

                    out.add(new MwvdExporter.ArcPiece(
                            arc,
                            pp.site1.id, pp.site2.id));
                }
            }
        }

        return out;
    }

    private MwvdModel.TrajectorySection sectionAt(MwvdModel.Trajectory traj, double t) {
        for (MwvdModel.TrajectorySection s : traj.secs) {
            if (s.includes(t)) return s;
        }
        return traj.secs.get(traj.secs.size() - 1);
    }
}