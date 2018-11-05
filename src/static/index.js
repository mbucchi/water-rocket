const bottle_volume = document.getElementById("bottleVol");
const fill_factor = document.getElementById("fillFactor");
const bottle_pressure = document.getElementById("bottlePressure");
const rocket_mass = document.getElementById("rocketMass");
const nozzle_area = document.getElementById("nozzleArea");
const launch_angle = document.getElementById("launchAngle");
const face_area = document.getElementById("faceArea");

document.getElementById("submit-form").addEventListener("submit", ev => {
  ev.preventDefault();

  fetch("/", {
    method: "POST",
    body: JSON.stringify({
      bottle_volume: Number.parseFloat(bottle_volume.value) * 1e-3,
      fill_factor: Number.parseFloat(fill_factor.value) * 1e-2,
      bottle_pressure: Number.parseFloat(bottle_pressure.value) * 1e5,
      rocket_mass: Number.parseFloat(rocket_mass.value) * 1e-3,
      nozzle_area: Number.parseFloat(nozzle_area.value) * 1e-4,
      launch_angle: Number.parseFloat(launch_angle.value),
      face_area: Number.parseFloat(face_area.value) * 1e-4
    })
  }).then(resp => {
    if (resp.status == 200) displayCharts();
    else alert("internal error running simulation");
  });
});

function displayCharts() {
  const rnd = new Date();
  document.getElementById(
    "trajectory-container"
  ).innerHTML = `<img width="100%" id="trajectory-img" src="/trajectory.png?${rnd}" alt="Trajectory">`;
  document.getElementById(
    "drag-container"
  ).innerHTML = `<img width="100%" id="drag-img" src="/drag.png?${rnd}" alt="Trajectory">`;
  document.getElementById(
    "thrust-container"
  ).innerHTML = `<img width="100%" id="thrust-img" src="/thrust.png?${rnd}" alt="Trajectory">`;
  document.getElementById(
    "water-speed-container"
  ).innerHTML = `<img width="100%" id="water-speed-img" src="/water-speed.png?${rnd}" alt="Trajectory">`;
  document.getElementById(
    "water-vol-container"
  ).innerHTML = `<img width="100%" id="water-vol-img" src="/water-vol.png?${rnd}" alt="Trajectory">`;

  fetch(`/animation?${rnd}`)
    .then(res => res.json())
    .then(displayAnimation);
}

function displayAnimation({ times, points }) {
  google.charts.load("current", { packages: ["corechart"] });
  google.charts.setOnLoadCallback(drawChart);

  function drawChart() {
    let maxValue =
      points.reduce((prev, [x, y]) => {
        let max = x > y ? x : y;
        return max > prev ? max : prev;
      }, 0) * 1.1;

    console.log(maxValue);

    let currentIdx = 0;
    let data;

    var options = {
      title: "Rocket Trajectory",
      hAxis: { title: "Distance (metres)", minValue: 0, maxValue },
      vAxis: { title: "Height (metres)", minValue: -1, maxValue },
      legend: "none"
    };

    var chart = new google.visualization.ScatterChart(
      document.getElementById("animation-container")
    );

    function setData() {
      data = google.visualization.arrayToDataTable([
        ["Height", "Distance"],
        points[currentIdx]
      ]);
      currentIdx += 3;
      drawChart();
      if (currentIdx >= points.length) {
        currentIdx = 0;
        setTimeout(setData, 5000);
      } else {
        setTimeout(setData, (times[currentIdx] - times[currentIdx - 3]) * 1000);
      }
    }

    function drawChart() {
      chart.draw(data, options);
    }

    setData();
  }
}
