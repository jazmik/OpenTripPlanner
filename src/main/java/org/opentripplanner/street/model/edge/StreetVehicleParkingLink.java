package org.opentripplanner.street.model.edge;

import com.google.common.collect.Sets;
import java.util.Set;
import org.locationtech.jts.geom.LineString;
import org.opentripplanner.framework.i18n.I18NString;
import org.opentripplanner.framework.tostring.ToStringBuilder;
import org.opentripplanner.routing.api.request.request.VehicleParkingRequest;
import org.opentripplanner.routing.vehicle_parking.VehicleParking;
import org.opentripplanner.street.model.vertex.StreetVertex;
import org.opentripplanner.street.model.vertex.VehicleParkingEntranceVertex;
import org.opentripplanner.street.search.TraverseMode;
import org.opentripplanner.street.search.state.State;
import org.opentripplanner.street.search.state.StateEditor;

/**
 * This represents the connection between a street vertex and a vehicle parking vertex.
 */
public class StreetVehicleParkingLink extends Edge {

  private final VehicleParkingEntranceVertex vehicleParkingEntranceVertex;

  public StreetVehicleParkingLink(StreetVertex fromv, VehicleParkingEntranceVertex tov) {
    super(fromv, tov);
    vehicleParkingEntranceVertex = tov;
  }

  public StreetVehicleParkingLink(VehicleParkingEntranceVertex fromv, StreetVertex tov) {
    super(fromv, tov);
    vehicleParkingEntranceVertex = fromv;
  }

  @Override
  public String toString() {
    return ToStringBuilder.of(this.getClass()).addObj("fromv", fromv).addObj("tov", tov).toString();
  }

  public State traverse(State s0) {
    // Disallow traversing two StreetBikeParkLinks in a row.
    // Prevents router using bike rental stations as shortcuts to get around
    // turn restrictions.
    if (s0.getBackEdge() instanceof StreetVehicleParkingLink) {
      return null;
    }

    var entrance = vehicleParkingEntranceVertex.getParkingEntrance();
    if (s0.getNonTransitMode() == TraverseMode.CAR) {
      if (!entrance.isCarAccessible()) {
        return null;
      }
    } else if (!entrance.isWalkAccessible()) {
      return null;
    }

    var vehicleParking = vehicleParkingEntranceVertex.getVehicleParking();
    final VehicleParkingRequest parkingRequest = s0.getRequest().parking();
    if (
      hasMissingRequiredTags(vehicleParking, parkingRequest.requiredTags()) ||
      hasBannedTags(vehicleParking, parkingRequest.bannedTags())
    ) {
      return null;
    }

    StateEditor s1 = s0.edit(this);

    if (isUnpreferredParking(parkingRequest, vehicleParking)) {
      s1.incrementWeight(parkingRequest.unpreferredTagCost());
    }

    s1.incrementWeight(1);
    s1.setBackMode(null);
    return s1.makeState();
  }

  @Override
  public I18NString getName() {
    return vehicleParkingEntranceVertex.getName();
  }

  public LineString getGeometry() {
    return null;
  }

  public double getDistanceMeters() {
    return 0;
  }

  private boolean hasBannedTags(VehicleParking vehicleParking, Set<String> bannedTags) {
    if (bannedTags.isEmpty()) {
      return false;
    }

    return vehicleParking.getTags().stream().anyMatch(bannedTags::contains);
  }

  private boolean hasMissingRequiredTags(VehicleParking vehicleParking, Set<String> requiredTags) {
    if (requiredTags.isEmpty()) {
      return false;
    }
    return !vehicleParking.getTags().containsAll(requiredTags);
  }

  private static boolean isUnpreferredParking(VehicleParkingRequest req, VehicleParking parking) {
    final var preferredTags = req.preferredTags();
    if (preferredTags.isEmpty()) {
      return false;
    } else {
      return Sets.intersection(preferredTags, parking.getTags()).isEmpty();
    }
  }
}
