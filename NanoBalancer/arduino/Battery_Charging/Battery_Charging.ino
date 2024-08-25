
#include <ArduinoMotorCarrier.h>
#include <Arduino_PMIC.h>
//Variable to store the battery voltage
float batteryVoltage;
int led = 13;

// the setup function runs once when you press reset or power the board
void setup() {
  Serial.begin(115200);
  pinMode(led, OUTPUT);
  
  if (controller.begin())
  {
    Serial.print("Nano Motor Shield connected, firmware version ");
    Serial.println(controller.getFWVersion());
  }
  else
  {
    Serial.println("Couldn't connect! Is the red led blinking? You may need to update the firmware with FWUpdater sketch");
    while (1);
  }
  // init battery settings
  if (!PMIC.begin()) {
    Serial.println("Failed to initialize PMIC!");
    while (1);
  }
  // pmic settings
  PMIC.setInputCurrentLimit(2.0);
  PMIC.setInputVoltageLimit(3.7);
  PMIC.setMinimumSystemVoltage(3.3);
  PMIC.setChargeVoltage(4.11);
  PMIC.setChargeCurrent(0.5);
  PMIC.setPreChargeCurrent(0.5);
  PMIC.setTermChargeCurrent(0.25);
  PMIC.disableWatchdog();
  
  // Enable the Charger
  PMIC.enableCharge();

  Serial.println("Initialization done!");
}

// the loop function runs over and over again forever
void loop() {
  digitalWrite(led, !digitalRead(led));    // Blink
  batteryVoltage = battery.getRaw()/236.0;
  Serial.print("Battery voltage: ");
  Serial.print(batteryVoltage,3);
  Serial.println("V");
  delay(1000); //wait for a few seconds
  
}
