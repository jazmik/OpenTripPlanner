package org.opentripplanner.apis.gtfs.model;

import org.locationtech.jts.geom.Coordinate;

/**
 * Class that contains a distance and geometry for each stop along a pattern
 */
public class PatternHopGeometryModel {

  private final  Double distance;

  private final Coordinate[] coordinates;

  public PatternHopGeometryModel(Double distance, Coordinate[] coordinates) {
    this.distance = distance;
    this.coordinates = coordinates;
  }

  public Coordinate[] getCoordinates() {
    return coordinates;
  }

  public Double getDistance() {
    return distance;
  }
}
